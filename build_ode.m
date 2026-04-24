function [ode_fun, out_fun] = build_ode(p)
%BUILD_ODE Build CasADi ODE and output map for raceway model.
%   [ode_fun, out_fun] = build_ode(p)
%   - ode_fun(x,u,d) returns xdot for x=[Xalg XO2 DIC Cat H T V].'
%   - out_fun(x,u,d) returns y=[pH DO Xalg_gL Depth].'
%
% Inputs:
%   p : parameter struct. Missing fields are filled with defaults.

import casadi.*

% if nargin < 1 || isempty(p)
%     p = struct();
% end
% p = create_p();

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
CO2_aq = DIC * H^2 / (Delta + p.reg_eps);
HCO3 = DIC * H * K1 / (Delta + p.reg_eps);
CO3 = DIC * K1 * K2 / (Delta + p.reg_eps);

% Electrically coupled proton dynamics terms
g_ratio = HCO3 / (DIC + p.reg_eps);
h_ratio = CO3  / (DIC + p.reg_eps);
dDelta_dH = 2*H + K1;
dg_dH = K1*(K1*K2 - H^2) / (Delta + p.reg_eps)^2;
dh_dH = -K1*K2 * dDelta_dH / (Delta + p.reg_eps)^2;

fDIC = -(g_ratio + 2*h_ratio);
fCat = 1;
fH = 1 - DIC*(dg_dH + 2*dh_dH) + KW/(H + p.reg_eps)^2;

% Gas transfer
Ug_CO2 = QCO2 / p.Asump; % equation (34)
Ug_O2  = Qair  / p.Asump; % equation (35)

kLa_CO2 = p.alpha_CO2 * (Ug_CO2)^p.beta_CO2; % equation (34)
kLa_O2  = p.alpha_O2  * (Ug_O2)^p.beta_O2; % equation (35)

kLa_CO2_eff = kLa_CO2 * p.Asump / (V + p.reg_eps); % equation (36)
kLa_O2_eff  = kLa_O2  * p.Asump / (V + p.reg_eps); % equation (37)

k_strip_CO2_eff = p.k_strip_CO2_O2 * kLa_O2_eff; % equation (38)
k_strip_O2_eff  = p.k_strip_O2_CO2 * kLa_CO2_eff; % equation (39)

Depth_val = (V - p.Vsump) / (p.W * p.L); % equation (6)

% Biological model
PAR = 0.46 * 4.56 * RadG; % equation (2)
Iav = PAR / (p.Ka * Depth_val * Xalg + p.reg_eps) * (1 - exp(-p.Ka * Depth_val * Xalg)); % equation (24)
mu_I = Iav^p.n_hill / (p.Ik^p.n_hill + Iav^p.n_hill); % equation (25)

mu_T = i_smooth_window(T, p.T_min, p.T_opt, p.T_max, p.reg_eps);
pH_sym = -log(H / 1000) / log(10);
mu_pH = i_smooth_window(pH_sym, p.pH_min, p.pH_opt, p.pH_max, p.reg_eps);

DO_sat = 100 * XO2 / (Xeq_O2 + p.reg_eps); % equation (27)
mu_DO = fmax(0, 1 - (DO_sat / p.DO_max)^p.m_DO); % equation (28)

P_rate = p.mu_max * mu_I * mu_T * mu_pH * mu_DO; % equation (29)
mu_g = p.eta_X * P_rate; % equation (30)
m_resp = p.m_min * (1 + p.k_resp_I * (1 - mu_I)) * p.Q10^((T - 20) / 10); % equation (31)

% Thermal model
Q_irrad = p.alpha_rad * p.A * RadG; % equation (9)

e_sat = 611.2 * exp(17.67 * T / (T + 243.5));
e_air = (RH/100) * 611.2 * exp(17.67 * Text / (Text + 243.5));
Cs_vap = 0.622 * e_sat / ((p.p_atm * 101325) - e_sat + p.reg_eps) * p.rho_w;
Ca_vap = 0.622 * e_air / ((p.p_atm * 101325) - e_air + p.reg_eps) * p.rho_w;

m_dot_e = p.k_m * p.A * fmax(Cs_vap - Ca_vap, 0); % equation (13)
lv = (2.501 - 0.00237 * T) * 1e6;
Q_evap = -lv * m_dot_e; % equation (14)
V_dot_e = m_dot_e / p.rho_w;

eps_sky = 0.70 + 5.95e-5 * (RH/100 * 4.596 * exp(17.27*Text/(237.3+Text)));
T_sky_K = (eps_sky * (Text + 273.15)^4)^(1/4); % ???????
Q_rad = p.sigma * p.eps_w * p.A * (T_sky_K^4 - T_K^4); % equation (10)

R_pp = p.x_liner / p.K_liner + p.x_ug / p.K_ug; % equation (11)
Q_cond = (p.A / (R_pp + p.reg_eps)) * (p.T_g - T); % equation (12)
Q_conv = p.h_c * p.A * (Text - T); % equation (15)

Q_dil = p.rho_w * p.Cp_w * Qd * p.T_in_medium; % equation (16)
Q_harv = -p.rho_w * p.Cp_w * Qh * T; % equation (17)
Q_mix = p.P_mix; % equation (18)

Cw_HX = p.rho_w * p.Cp_w * Qw; % equation (21)

eps_HX = 1 - exp(-p.UA / (Cw_HX + p.reg_eps)); % equation (22)
Q_HX = Cw_HX * (Tin_hx - T) * eps_HX; % equation (23)

Q_Sigma = Q_irrad + Q_rad + Q_cond + Q_evap + Q_conv + Q_dil + Q_harv + Q_mix + Q_HX; % equation (8)

% Balances
V_dot = Qd - Qh - V_dot_e; % equation (32)

T_dot = Q_Sigma / (p.rho_w * p.Cp_w * V + p.reg_eps) - (T / (V + p.reg_eps)) * V_dot;
Xalg_dot = (mu_g - m_resp) * Xalg - (Qd / (V + p.reg_eps)) * Xalg;

DIC_dot = (Qd / (V + p.reg_eps)) * (p.DIC_in - DIC) ...
        - P_rate * Xalg / (p.Y_CO2 / p.M_CO2) ...
        + m_resp * Xalg / (p.Y_CO2 / p.M_CO2) ...
        + kLa_CO2_eff * (COinj_CO2 - CO2_aq) ...
        + p.k_atm_CO2 * (COeq_CO2 - CO2_aq) ...
        + p.k_pw_CO2 * p.W * p.Lpw * Depth_val / (V + p.reg_eps) * (COeq_CO2 - CO2_aq) ...
        + k_strip_CO2_eff * (COeq_CO2 - CO2_aq);

Cat_dot = (Qd / (V + p.reg_eps)) * (p.Cat_in - Cat);

XO2_dot = (Qd / (V + p.reg_eps)) * (Xeq_O2 - XO2) ...
        + P_rate * Xalg / (p.Y_O2 / p.M_O2) ...
        - m_resp * Xalg / (p.Y_O2 / p.M_O2) ...
        + kLa_O2_eff * (Xeq_O2 - XO2) ...
        + p.k_atm_O2 * (Xeq_O2 - XO2) ...
        + p.k_pw_O2 * p.W * p.Lpw * Depth_val / (V + p.reg_eps) * (Xeq_O2 - XO2) ...
        - k_strip_O2_eff * XO2;

H_dot = -(fDIC * DIC_dot + fCat * Cat_dot) / (fH + p.reg_eps);

xdot = vertcat(Xalg_dot, XO2_dot, DIC_dot, Cat_dot, H_dot, T_dot, V_dot);

% Output map
pH_out = -log(H / 1000) / log(10);
DO_out = 100 * XO2 / (Xeq_O2 + p.reg_eps);
Xalg_gL = Xalg / 1000;
Depth_out = (V - p.Vsump) / (p.W * p.L + p.reg_eps);

y = vertcat(pH_out, DO_out, Xalg_gL, Depth_out);

ode_fun = Function('ode_fun', {x, u, d}, {xdot}, {'x','u','d'}, {'xdot'});
out_fun = Function('out_fun', {x, u, d}, {y}, {'x','u','d'}, {'y'});

end

function mu = i_smooth_window(val, a, b, c, reg_eps)
import casadi.*

r_left = (val - a) / (b - a + reg_eps);
mu_left = 3*r_left^2 - 2*r_left^3;

r_right = (c - val) / (c - b + reg_eps);
mu_right = 3*r_right^2 - 2*r_right^3;

mu = if_else(val <= a, 0, if_else(val <= b, mu_left, if_else(val < c, mu_right, 0)));
mu = fmin(fmax(mu, 0), 1);

end


