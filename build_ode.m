function [ode_fun, out_fun] = build_ode(p)
%BUILD_ODE Build CasADi ODE and output map for raceway model.
%   [ode_fun, out_fun] = build_ode(p)
%   - ode_fun(x,u,d) returns xdot for x=[Xalg XO2 DIC Cat H T V].'
%   - out_fun(x,u,d) returns y=[pH DO Xalg_gL Depth].'
%
% Inputs:
%   p : parameter struct. Missing fields are filled with defaults.

import casadi.*

if nargin < 1 || isempty(p)
    p = struct();
end
p = i_fill_defaults(p);

% Symbolic variables
x = MX.sym('x', 7);   % [Xalg, XO2, DIC, Cat, H, T, V]
u = MX.sym('u', 6);   % [QCO2, Qair, Qd, Qh, Qw, Tin_hx]
d = MX.sym('d', 4);   % [RadG, RH, Text, Uwind]

% States
Xalg = x(1);
XO2  = x(2);
DIC  = x(3);
Cat  = x(4);
H    = x(5);
T    = x(6);
V    = x(7);

% Inputs
QCO2   = u(1);
Qair   = u(2);
Qd     = u(3);
Qh     = u(4);
Qw     = u(5);
Tin_hx = u(6);

% Disturbances
RadG = d(1);
RH   = d(2);
Text = d(3);

T_K = T + 273.15;

% Gas-liquid equilibria
KH_O2  = p.KH_O2_ref  * exp(p.C_O2  * (1/T_K - 1/p.T_ref));
KH_CO2 = p.KH_CO2_ref * exp(p.C_CO2 * (1/T_K - 1/p.T_ref));

Xeq_O2   = KH_O2 * p.p_atm * p.y_O2;
COeq_CO2 = KH_CO2 * p.p_atm * p.y_CO2_atm;
COinj_CO2 = KH_CO2 * p.p_atm * p.y_CO2_inj;

K1 = p.K1_ref * exp(-p.dH_K1/p.R_gas * (1/T_K - 1/p.T_ref));
K2 = p.K2_ref * exp(-p.dH_K2/p.R_gas * (1/T_K - 1/p.T_ref));
KW = p.KW_ref * exp(-p.dH_KW/p.R_gas * (1/T_K - 1/p.T_ref));

% Carbonate species
Delta = H^2 + H*K1 + K1*K2;
CO2_aq = DIC * H^2 / (Delta + 1e-30);
HCO3 = DIC * H * K1 / (Delta + 1e-30);
CO3 = DIC * K1 * K2 / (Delta + 1e-30);

% Electrically coupled proton dynamics terms
g_ratio = HCO3 / (DIC + 1e-30);
h_ratio = CO3  / (DIC + 1e-30);
dDelta_dH = 2*H + K1;
dg_dH = K1*(K1*K2 - H^2) / (Delta + 1e-30)^2;
dh_dH = -K1*K2 * dDelta_dH / (Delta + 1e-30)^2;

fDIC = -(g_ratio + 2*h_ratio);
fCat = 1;
fH = 1 - DIC*(dg_dH + 2*dh_dH) + KW/(H + 1e-30)^2;

% Gas transfer
Ug_CO2 = QCO2 / (p.Asump + 1e-30);
Ug_O2  = Qair  / (p.Asump + 1e-30);

kLa_CO2 = p.alpha_CO2 * (Ug_CO2 + 1e-10)^p.beta_CO2;
kLa_O2  = p.alpha_O2  * (Ug_O2  + 1e-10)^p.beta_O2;

kLa_CO2_eff = kLa_CO2 * p.Asump / (V + 1e-30);
kLa_O2_eff  = kLa_O2  * p.Asump / (V + 1e-30);

k_strip_CO2_eff = p.k_strip_CO2_O2 * kLa_O2_eff;
k_strip_O2_eff  = p.k_strip_O2_CO2 * kLa_CO2_eff;

Depth_val = (V - p.Vsump) / (p.W * p.L + 1e-30);

% Biological model
PAR = 0.46 * 4.56 * RadG;
Iav = PAR / (p.Ka * Depth_val * Xalg + 1e-30) * (1 - exp(-p.Ka * Depth_val * Xalg));
mu_I = Iav^p.n_hill / (p.Ik^p.n_hill + Iav^p.n_hill + 1e-30);

mu_T = i_smooth_window(T, p.T_min, p.T_opt, p.T_max);
pH_sym = -log(H / 1000) / log(10);
mu_pH = i_smooth_window(pH_sym, p.pH_min, p.pH_opt, p.pH_max);

DO_sat = 100 * XO2 / (Xeq_O2 + 1e-30);
mu_DO = fmax(0, 1 - (DO_sat / p.DO_max)^p.m_DO);

P_rate = p.mu_max * mu_I * mu_T * mu_pH * mu_DO;
mu_g = p.eta_X * P_rate;
m_resp = p.m_min * (1 + p.k_resp_I * (1 - mu_I)) * p.Q10^((T - 20) / 10);

% Thermal model
Q_irrad = p.alpha_rad * p.A * RadG;

e_sat = 611.2 * exp(17.67 * T / (T + 243.5));
e_air = (RH/100) * 611.2 * exp(17.67 * Text / (Text + 243.5));
Cs_vap = 0.622 * e_sat / ((p.p_atm * 101325) - e_sat + 1e-30) * p.rho_w;
Ca_vap = 0.622 * e_air / ((p.p_atm * 101325) - e_air + 1e-30) * p.rho_w;

m_dot_e = p.k_m * p.A * fmax(Cs_vap - Ca_vap, 0);
lv = (2.501 - 0.00237 * T) * 1e6;
Q_evap = -lv * m_dot_e;
V_dot_e = m_dot_e / p.rho_w;

eps_sky = 0.70 + 5.95e-5 * (RH/100 * 4.596 * exp(17.27*Text/(237.3+Text)));
T_sky_K = (eps_sky * (Text + 273.15)^4)^(1/4);
Q_rad = p.sigma * p.eps_w * p.A * (T_sky_K^4 - T_K^4);

R_pp = p.x_liner / p.K_liner + p.x_ug / p.K_ug;
Q_cond = (p.A / (R_pp + 1e-30)) * (p.T_g - T);
Q_conv = p.h_c * p.A * (Text - T);

Q_dil = p.rho_w * p.Cp_w * Qd * p.T_in_medium;
Q_harv = -p.rho_w * p.Cp_w * Qh * T;
Q_mix = p.P_mix;

Cw_HX = p.rho_w * p.Cp_w * Qw;
eps_HX = 1 - exp(-p.UA_HX / (Cw_HX + 1e-30));
Q_HX = Cw_HX * (Tin_hx - T) * eps_HX;

Q_Sigma = Q_irrad + Q_rad + Q_cond + Q_evap + Q_conv + Q_dil + Q_harv + Q_mix + Q_HX;

% Balances
V_dot = Qd - Qh - V_dot_e;

T_dot = Q_Sigma / (p.rho_w * p.Cp_w * V + 1e-30) - (T / (V + 1e-30)) * V_dot;
Xalg_dot = (mu_g - m_resp) * Xalg - (Qd / (V + 1e-30)) * Xalg;

DIC_dot = (Qd / (V + 1e-30)) * (p.DIC_in - DIC) ...
        - P_rate * Xalg / (p.Y_CO2 / p.M_CO2) ...
        + m_resp * Xalg / (p.Y_CO2 / p.M_CO2) ...
        + kLa_CO2_eff * (COinj_CO2 - CO2_aq) ...
        + p.k_atm_CO2 * (COeq_CO2 - CO2_aq) ...
        + p.k_pw_CO2 * p.W * p.Lpw * Depth_val / (V + 1e-30) * (COeq_CO2 - CO2_aq) ...
        + k_strip_CO2_eff * (COeq_CO2 - CO2_aq);

Cat_dot = (Qd / (V + 1e-30)) * (p.Cat_in - Cat);

XO2_dot = (Qd / (V + 1e-30)) * (Xeq_O2 - XO2) ...
        + P_rate * Xalg / (p.Y_O2 / p.M_O2) ...
        - m_resp * Xalg / (p.Y_O2 / p.M_O2) ...
        + kLa_O2_eff * (Xeq_O2 - XO2) ...
        + p.k_atm_O2 * (Xeq_O2 - XO2) ...
        + p.k_pw_O2 * p.W * p.Lpw * Depth_val / (V + 1e-30) * (Xeq_O2 - XO2) ...
        - k_strip_O2_eff * XO2;

H_dot = -(fDIC * DIC_dot + fCat * Cat_dot) / (fH + 1e-30);

xdot = vertcat(Xalg_dot, XO2_dot, DIC_dot, Cat_dot, H_dot, T_dot, V_dot);

% Output map
pH_out = -log(H / 1000) / log(10);
DO_out = 100 * XO2 / (Xeq_O2 + 1e-30);
Xalg_gL = Xalg / 1000;
Depth_out = (V - p.Vsump) / (p.W * p.L + 1e-30);

y = vertcat(pH_out, DO_out, Xalg_gL, Depth_out);

ode_fun = Function('ode_fun', {x, u, d}, {xdot}, {'x','u','d'}, {'xdot'});
out_fun = Function('out_fun', {x, u, d}, {y}, {'x','u','d'}, {'y'});

end

function mu = i_smooth_window(val, a, b, c)
import casadi.*

r_left = (val - a) / (b - a + 1e-30);
mu_left = 3*r_left^2 - 2*r_left^3;

r_right = (c - val) / (c - b + 1e-30);
mu_right = 3*r_right^2 - 2*r_right^3;

mu = if_else(val <= a, 0, if_else(val <= b, mu_left, if_else(val < c, mu_right, 0)));
mu = fmin(fmax(mu, 0), 1);

end

function p = i_fill_defaults(p)
% Geometry
if ~isfield(p,'L'), p.L = 80; end                              % [m]         channel length
if ~isfield(p,'W'), p.W = 1; end                               % [m]         channel width
if ~isfield(p,'A'), p.A = p.L * p.W; end                      % [m²]        illuminated surface
if ~isfield(p,'rsump'), p.rsump = 0.4; end                     % [m]         sump radius
if ~isfield(p,'hsump'), p.hsump = 1.12; end                    % [m]         sump depth
if ~isfield(p,'Asump'), p.Asump = pi * p.rsump^2; end          % [m²]        sump cross-section
if ~isfield(p,'Vsump'), p.Vsump = p.Asump * p.hsump; end       % [m³]        sump volume

% Water and thermal
if ~isfield(p,'rho_w'), p.rho_w = 997; end                    % [kg/m³]     water density
if ~isfield(p,'Cp_w'), p.Cp_w = 4182; end                      % [J/(kg·K)]  water heat capacity
if ~isfield(p,'alpha_rad'), p.alpha_rad = 0.7; end             % [-]         shortwave absorptance
if ~isfield(p,'sigma'), p.sigma = 5.6697e-8; end               % [W/(m²·K⁴)] Stefan-Boltzmann constant
if ~isfield(p,'eps_w'), p.eps_w = 0.97; end                    % [-]         water emissivity
if ~isfield(p,'x_liner'), p.x_liner = 0.002; end               % [m]         liner thickness
if ~isfield(p,'K_liner'), p.K_liner = 0.33; end                % [W/(m·K)]   liner conductivity
if ~isfield(p,'x_ug'), p.x_ug = 0.20; end                      % [m]         subgrade thickness
if ~isfield(p,'K_ug'), p.K_ug = 1.0; end                       % [W/(m·K)]   subgrade conductivity
if ~isfield(p,'T_g'), p.T_g = 17.4; end                        % [°C]        ground reference temperature
if ~isfield(p,'k_m'), p.k_m = 0.01; end                        % [m/s]       evaporation mass-transfer coeff.
if ~isfield(p,'h_c'), p.h_c = 10; end                          % [W/(m²·K)]  convective heat-transfer coeff.
if ~isfield(p,'P_mix'), p.P_mix = 200; end                      % [W]         paddlewheel mechanical power
if ~isfield(p,'UA_HX'), p.UA_HX = 500; end                     % [W/K]       heat exchanger overall UA
if ~isfield(p,'T_in_medium'), p.T_in_medium = 20; end          % [°C]        dilution medium temperature
if ~isfield(p,'Lpw'), p.Lpw = 2; end                           % [m]         paddlewheel zone effective length

% Biology
if ~isfield(p,'mu_max'), p.mu_max = 0.8/86400; end             % [s⁻¹]       max. specific photosynthesis (0.8 d⁻¹)
if ~isfield(p,'Ka'), p.Ka = 0.2; end                           % [m²/g]      light extinction coefficient
if ~isfield(p,'Ik'), p.Ik = 200; end                           % [µmol/(m²·s)] half-sat. irradiance
if ~isfield(p,'n_hill'), p.n_hill = 1; end                      % [-]         Hill exponent
if ~isfield(p,'T_min'), p.T_min = 5; end                        % [°C]        temperature window minimum
if ~isfield(p,'T_opt'), p.T_opt = 30; end                      % [°C]        optimal temperature
if ~isfield(p,'T_max'), p.T_max = 40; end                      % [°C]        temperature window maximum
if ~isfield(p,'pH_min'), p.pH_min = 6; end                      % [-]         pH window minimum
if ~isfield(p,'pH_opt'), p.pH_opt = 8; end                     % [-]         optimal pH
if ~isfield(p,'pH_max'), p.pH_max = 10; end                    % [-]         pH window maximum
if ~isfield(p,'DO_max'), p.DO_max = 383.21; end                % [%]         DO inhibition threshold
if ~isfield(p,'m_DO'), p.m_DO = 2; end                         % [-]         inhibition exponent
if ~isfield(p,'eta_X'), p.eta_X = 0.6; end                     % [-]         gross-to-growth allocation
if ~isfield(p,'Y_CO2'), p.Y_CO2 = 44/30; end                 % [gCO2/gXalg] CO2 yield coefficient
if ~isfield(p,'Y_O2'), p.Y_O2 = 32/30; end                   % [gO2/gXalg]  O2 yield coefficient
if ~isfield(p,'M_CO2'), p.M_CO2 = 44; end                      % [g/mol]     molar mass of CO2
if ~isfield(p,'M_O2'), p.M_O2 = 32; end                       % [g/mol]     molar mass of O2
if ~isfield(p,'m_min'), p.m_min = 0.02/86400; end              % [s⁻¹]       basal respiration rate at 20°C
if ~isfield(p,'k_resp_I'), p.k_resp_I = 0.5; end               % [-]         low-light maintenance gain
if ~isfield(p,'Q10'), p.Q10 = 2.0; end                         % [-]         temperature sensitivity (Q10)

% Chemistry and transfer
if ~isfield(p,'T_ref'), p.T_ref = 298.15; end                  % [K]         reference temperature
if ~isfield(p,'R_gas'), p.R_gas = 8.314; end                   % [J/(mol·K)] ideal gas constant
if ~isfield(p,'p_atm'), p.p_atm = 1.0; end                     % [atm]       atmospheric pressure
if ~isfield(p,'y_O2'), p.y_O2 = 0.2095; end                      % [-]         O2 molar fraction in air
if ~isfield(p,'y_CO2_atm'), p.y_CO2_atm = 0.00039; end         % [-]         CO2 molar fraction in air (~380 ppm)
if ~isfield(p,'y_CO2_inj'), p.y_CO2_inj = 1.0; end            % [-]         CO2 molar fraction in injection gas (pure CO2)
if ~isfield(p,'KH_O2_ref'), p.KH_O2_ref = 1.3e-3 * 1e3; end   % [mol/(m³·atm)] Henry const. for O2 at T_ref
if ~isfield(p,'KH_CO2_ref'), p.KH_CO2_ref = 3.4e-2 * 1e3; end % [mol/(m³·atm)] Henry const. for CO2 at T_ref
if ~isfield(p,'C_O2'), p.C_O2 = 1700; end                      % [K]         Van't Hoff constant for O2
if ~isfield(p,'C_CO2'), p.C_CO2 = 2400; end                   % [K]         Van't Hoff constant for CO2
if ~isfield(p,'K1_ref'), p.K1_ref = 10^(-6.3) * 1e3; end         % [mol/m³]    first carbonate dissociation const.
if ~isfield(p,'K2_ref'), p.K2_ref = 10^(-10.33) * 1e3; end        % [mol/m³]    second carbonate dissociation const.
if ~isfield(p,'KW_ref'), p.KW_ref = 1e-14 * 1e6; end          % [mol²/m⁶]   water autoprotolysis const.
if ~isfield(p,'dH_K1'), p.dH_K1 = -7500; end                   % [J/mol]     enthalpy of K1
if ~isfield(p,'dH_K2'), p.dH_K2 = -15000; end                 % [J/mol]     enthalpy of K2
if ~isfield(p,'dH_KW'), p.dH_KW = 55900; end                  % [J/mol]     enthalpy of KW
if ~isfield(p,'alpha_CO2'), p.alpha_CO2 = 0.05; end            % [-]         kLa scale factor for CO2
if ~isfield(p,'beta_CO2'), p.beta_CO2 = 0.7; end               % [-]         kLa power-law exponent for CO2
if ~isfield(p,'alpha_O2'), p.alpha_O2 = 0.04; end              % [-]         kLa scale factor for O2
if ~isfield(p,'beta_O2'), p.beta_O2 = 0.7; end                 % [-]         kLa power-law exponent for O2
if ~isfield(p,'k_strip_CO2_O2'), p.k_strip_CO2_O2 = 0.05; end  % [-]         CO2 cross-stripping factor
if ~isfield(p,'k_strip_O2_CO2'), p.k_strip_O2_CO2 = 0.02; end  % [-]         O2 cross-stripping factor
if ~isfield(p,'k_atm_CO2'), p.k_atm_CO2 = 1e-5; end            % [s⁻¹]       atmospheric CO2 exchange coeff.
if ~isfield(p,'k_atm_O2'), p.k_atm_O2 = 2e-5; end              % [s⁻¹]       atmospheric O2 exchange coeff.
if ~isfield(p,'k_pw_CO2'), p.k_pw_CO2 = 1e-5; end              % [s⁻¹]       paddlewheel CO2 exchange coeff.
if ~isfield(p,'k_pw_O2'), p.k_pw_O2 = 2e-5; end                % [s⁻¹]       paddlewheel O2 exchange coeff.
if ~isfield(p,'DIC_in'), p.DIC_in = 2; end                     % [mol/m³]    inlet dissolved inorganic carbon
if ~isfield(p,'Cat_in'), p.Cat_in = 5; end                     % [mol/m³]    inlet strong cations

end
