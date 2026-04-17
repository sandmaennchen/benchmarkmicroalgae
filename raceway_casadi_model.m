%% =========================================================================
%  Microalgae Raceway Reactor – CasADi MATLAB Model
%  Based on: Rodríguez-Miranda et al. (2025), arXiv:2512.15916v1
%
%  Model equations implemented:
%    - Thermal submodel      (Section 2.3)
%    - Biological submodel   (Section 2.4)
%    - Engineering submodel  (Section 2.5)
%    - Gas-liquid equilibria (Section 2.6)
%
%  Usage:
%    Run this script directly for a demo open-loop simulation, or call
%    build_raceway_ode() from your own script/optimiser to obtain CasADi
%    Function objects for the ODE right-hand side and the output map.
%
%  Dependencies:
%    CasADi (https://web.casadi.org/) must be on the MATLAB path.
%
%  State vector  x = [Xalg, XO2, DIC, Cat, H, T, V]   (7 states)
%  Input vector  u = [QCO2, Qair, Qd, Qh, Qw, Tin_hx]  (6 inputs)
%  Disturbance   d = [RadG, RH, Text, Uwind]             (4 disturbances)
%  Output vector y = [pH, DO, Xalg_gL, Depth]           (4 outputs)
%
%  Author note: equation numbers in comments refer to the paper above.
% =========================================================================

import casadi.*

%% -------------------------------------------------------------------------
%  1.  PARAMETERS
%  All numerical values are taken from Table 1 and Section 2 of the paper.
%  Values not given explicitly are set to physically reasonable defaults
%  (marked with *default*) – replace with your calibrated values.
% -------------------------------------------------------------------------

p = struct();

% --- Reactor geometry ---
p.L         = 40;           % [m]   channel length (two 40-m channels)
p.W         = 1;            % [m]   channel width
p.A         = p.W * p.L;   % [m²]  illuminated surface = 80 m²
p.rsump     = 0.4;          % [m]   sump radius
p.hsump     = 1.12;         % [m]   sump depth
p.Asump     = pi * p.rsump^2;          % [m²]  sump cross-section
p.Vsump     = p.Asump * p.hsump;       % [m³]  sump volume

% --- Water physical properties ---
p.rho_w     = 1000;         % [kg/m³]
p.Cp_w      = 4182;         % [J/(kg·K)]

% --- Thermal submodel parameters ---
p.alpha_rad = 0.7;          % [-]   shortwave absorptance (*default*)
p.sigma     = 5.67e-8;      % [W/(m²·K⁴)] Stefan–Boltzmann
p.eps_w     = 0.97;         % [-]   water emissivity (*default*)
p.x_liner   = 0.003;        % [m]   liner thickness (3 mm)
p.K_liner   = 0.35;         % [W/(m·K)] liner conductivity (*default*)
p.x_ug      = 0.20;         % [m]   subgrade thickness (*default*)
p.K_ug      = 1.0;          % [W/(m·K)] subgrade conductivity (*default*)
p.T_g       = 20;           % [°C]  ground reference temperature (*default*)
p.k_m       = 0.01;         % [m/s] evaporation mass-transfer coeff. (*default*)
p.h_c       = 10;           % [W/(m²·K)] convective heat-transfer coeff. (*default*)
p.P_mix     = 200;          % [W]   paddlewheel mechanical power (*default*)

% Heat exchanger (spiral tube in sump)
p.UA_HX     = 500;          % [W/K] overall UA (*default*)

% --- Biological submodel parameters ---
p.mu_max    = 0.8 / 86400;  % [s⁻¹] max. specific photosynthesis (0.8 d⁻¹)
p.Ka        = 0.2;          % [m²/g] light extinction (*default*)
p.Ik        = 200;          % [µmol/(m²·s)] half-sat. irradiance (*default*)
p.n_hill    = 1;            % [-]   Hill exponent (*default*)

% Temperature window  [Tmin, Topt, Tmax] °C
p.T_min     = 5;
p.T_opt     = 30;           % optimal temperature (paper Section 2.1)
p.T_max     = 40;

% pH window  [pHmin, pHopt, pHmax]
p.pH_min    = 6;
p.pH_opt    = 8;            % optimal pH (paper Section 2.1)
p.pH_max    = 10;

% DO inhibition
p.DO_max    = 383.21;       % [%]   inhibition threshold (paper Section 2.1)
p.m_DO      = 2;            % [-]   inhibition exponent (*default*)

% Biomass allocation & stoichiometry
p.eta_X     = 0.6;          % [-]   gross-to-growth allocation (*default*)
p.Y_CO2     = 1.92;         % [gCO2/gXalg] (*default*)
p.Y_O2      = 1.41;         % [gO2/gXalg]  (*default*)
p.M_CO2     = 44;           % [g/mol]
p.M_O2      = 32;           % [g/mol]

% Maintenance / respiration
p.m_min     = 0.02 / 86400; % [s⁻¹] basal respiration at 20°C (*default*)
p.k_resp_I  = 0.5;          % [-]   low-light maintenance gain (*default*)
p.Q10       = 2.0;          % [-]   temperature sensitivity (*default*)

% --- Gas-liquid equilibria (Henry's law, Van't Hoff) ---
p.T_ref     = 298.15;       % [K]
p.R_gas     = 8.314;        % [J/(mol·K)]
p.p_atm     = 1;            % [atm]
p.y_O2      = 0.21;         % [-]   O2 molar fraction in air
p.y_CO2_atm = 3.8e-4;       % [-]   CO2 molar fraction in air (~380 ppm)
p.y_CO2_inj = 1.0;          % [-]   pure CO2 injection

% Henry constants at T_ref [mol/(m³·atm)]
p.KH_O2_ref  = 1.3e-3 * 1e3;   % convert mol/(L·atm) → mol/(m³·atm)
p.KH_CO2_ref = 3.4e-2 * 1e3;
p.C_O2       = 1700;            % [K]  Van't Hoff constant for O2 (*default*)
p.C_CO2      = 2400;            % [K]  Van't Hoff constant for CO2 (*default*)

% Carbonate equilibrium at T_ref (25 °C) [mol/L → mol/m³ where needed]
p.K1_ref    = 4.3e-7 * 1e3;    % [mol/m³]  first dissociation
p.K2_ref    = 4.7e-11 * 1e3;   % [mol/m³]  second dissociation
p.KW_ref    = 1e-14 * 1e6;     % [mol²/m⁶] water autoprotolysis
p.dH_K1     = -7500;           % [J/mol]   (*default*)
p.dH_K2     = -15000;          % [J/mol]   (*default*)
p.dH_KW     = 55900;           % [J/mol]   (*default*)

% Gas-transfer coefficients (sump, power-law)
p.alpha_CO2 = 0.05;     % [–]   kLa scale for CO2 (*default*)
p.beta_CO2  = 0.7;      % [–]   kLa exponent for CO2 (*default*)
p.alpha_O2  = 0.04;     % [–]   kLa scale for O2 (*default*)
p.beta_O2   = 0.7;      % [–]   kLa exponent for O2 (*default*)

% Cross-stripping factors
p.k_strip_CO2_O2 = 0.05;  % (*default*)
p.k_strip_O2_CO2 = 0.02;  % (*default*)

% Atmospheric and paddlewheel exchange [s⁻¹]
p.k_atm_CO2 = 1e-5;   % (*default*)
p.k_atm_O2  = 2e-5;   % (*default*)
p.k_pw_CO2  = 1e-5;   % (*default*)
p.k_pw_O2   = 2e-5;   % (*default*)

% Inlet medium composition
p.DIC_in    = 2;    % [mol/m³] inlet DIC (*default*)
p.Cat_in    = 5;    % [mol/m³] inlet strong cations (*default*)
p.T_in_medium = 20; % [°C]    dilution medium temperature (*default*)

% Liquid velocity for sump dispersion
p.Lpw       = 2;    % [m]     paddlewheel zone effective length (*default*)

%% -------------------------------------------------------------------------
%  2.  BUILD CASADI ODE FUNCTION
% -------------------------------------------------------------------------

[ode_fun, out_fun] = build_raceway_ode(p);

fprintf('CasADi ODE function built successfully.\n');
fprintf('  States  (x): Xalg, XO2, DIC, Cat, H, T, V\n');
fprintf('  Inputs  (u): QCO2, Qair, Qd, Qh, Qw, Tin_hx\n');
fprintf('  Disturb (d): RadG, RH, Text, Uwind\n');
fprintf('  Outputs (y): pH, DO, Xalg_gL, Depth\n\n');

%% -------------------------------------------------------------------------
%  3.  DEMO: OPEN-LOOP SIMULATION (Euler integration, 24-hour day)
% -------------------------------------------------------------------------

Tm     = 60;            % [s]  integration step (as in benchmark)
t_end  = 24 * 3600;     % [s]  simulate 24 hours
N_sim  = t_end / Tm;

% --- Initial conditions ---
x0 = [ 500;    % Xalg  [g/m³]  = 0.5 g/L
        0.26;  % XO2   [mol/m³] ≈ 100% saturation at 20°C
        2.0;   % DIC   [mol/m³]
        5.0;   % Cat   [mol/m³]
        1e-8;  % H     [mol/m³] → pH ≈ 8
        20;    % T     [°C]
        12 ];  % V     [m³]  ≈ depth 0.15 m × 80 m²

% --- Constant inputs (all in m³/s) ---
u_const = [ 5e-5;   % QCO2   [m³/s]
            1e-4;   % Qair   [m³/s]
            0;      % Qd     [m³/s]
            0;      % Qh     [m³/s]
            0;      % Qw     [m³/s]  (HX off)
            20  ];  % Tin_hx [°C]

% --- Synthetic sinusoidal irradiance disturbance ---
t_vec  = (0:N_sim-1) * Tm;
RadG_v = max(0, 800 * sin(pi * t_vec / t_end));  % peak 800 W/m²
RH_v   = 50 * ones(1, N_sim);
Text_v = 20 + 10 * sin(pi * t_vec / t_end);
Uwind_v= 3 * ones(1, N_sim);

% --- Storage ---
X_log  = zeros(7, N_sim+1);
Y_log  = zeros(4, N_sim);
X_log(:,1) = x0;

x_k = x0;
for k = 1:N_sim
    d_k = [RadG_v(k); RH_v(k); Text_v(k); Uwind_v(k)];

    % Evaluate ODE
    xdot_k = full(ode_fun(x_k, u_const, d_k));

    % Evaluate outputs
    y_k = full(out_fun(x_k, u_const, d_k));
    Y_log(:, k) = y_k;

    % Forward Euler step
    x_next = x_k + Tm * xdot_k;

    % Hard lower bounds to avoid numerical blow-up in demo
    x_next(1) = max(x_next(1), 1e-6);   % Xalg > 0
    x_next(2) = max(x_next(2), 1e-9);   % XO2  > 0
    x_next(3) = max(x_next(3), 1e-9);   % DIC  > 0
    x_next(5) = max(x_next(5), 1e-12);  % H    > 0
    x_next(7) = max(x_next(7), p.Vsump + 0.01); % V > Vsump

    X_log(:, k+1) = x_next;
    x_k = x_next;
end

% --- Plot results ---
t_h = t_vec / 3600;

figure('Name','Raceway Open-Loop Demo','NumberTitle','off','Position',[100 100 900 800]);

subplot(5,1,1);
plot(t_h, Y_log(1,:), 'b', 'LineWidth', 1.5); hold on;
yline(8, 'k--', 'pH setpoint');
ylabel('pH [-]'); title('pH'); grid on; xlim([0 24]);

subplot(5,1,2);
plot(t_h, Y_log(2,:), 'r', 'LineWidth', 1.5); hold on;
yline(150, 'k--', 'DO setpoint');
ylabel('DO [%sat]'); title('Dissolved Oxygen'); grid on; xlim([0 24]);

subplot(5,1,3);
plot(t_h, Y_log(3,:), 'g', 'LineWidth', 1.5);
ylabel('X_{alg} [g/L]'); title('Biomass Concentration'); grid on; xlim([0 24]);

subplot(5,1,4);
plot(t_h, X_log(6, 1:N_sim), 'm', 'LineWidth', 1.5); hold on;
yline(30, 'k--', 'T setpoint');
ylabel('T [°C]'); title('Raceway Temperature'); grid on; xlim([0 24]);

subplot(5,1,5);
plot(t_h, RadG_v, 'k', 'LineWidth', 1.5);
ylabel('RadG [W/m²]'); title('Solar Irradiance (disturbance)');
xlabel('Time [h]'); grid on; xlim([0 24]);

sgtitle('Microalgae Raceway – CasADi Open-Loop Demo (24 h)');

fprintf('Demo simulation complete. See figure for trajectories.\n');

%% =========================================================================
%  LOCAL FUNCTION: build_raceway_ode
%  Returns two CasADi Functions:
%    ode_fun(x, u, d)  →  xdot   (7 × 1)
%    out_fun(x, u, d)  →  y      (4 × 1)
% =========================================================================

function [ode_fun, out_fun] = build_raceway_ode(p)

import casadi.*

% ---- Symbolic variables --------------------------------------------------
x   = MX.sym('x', 7);   % [Xalg, XO2, DIC, Cat, H, T, V]
u   = MX.sym('u', 6);   % [QCO2, Qair, Qd, Qh, Qw, Tin_hx]
d   = MX.sym('d', 4);   % [RadG, RH, Text, Uwind]

% ---- Unpack states -------------------------------------------------------
Xalg = x(1);   % [g/m³]
XO2  = x(2);   % [mol/m³]
DIC  = x(3);   % [mol/m³]
Cat  = x(4);   % [mol/m³]
H    = x(5);   % [mol/m³]  (H⁺ concentration)
T    = x(6);   % [°C]
V    = x(7);   % [m³]

% ---- Unpack inputs -------------------------------------------------------
QCO2   = u(1);  % [m³/s]
Qair   = u(2);  % [m³/s]
Qd     = u(3);  % [m³/s]  dilution inflow
Qh     = u(4);  % [m³/s]  harvest outflow
Qw     = u(5);  % [m³/s]  HX water flow
Tin_hx = u(6);  % [°C]    HX inlet temperature

% ---- Unpack disturbances -------------------------------------------------
RadG   = d(1);  % [W/m²]
RH     = d(2);  % [%]
Text   = d(3);  % [°C]
% Uwind = d(4); % [m/s]  (not used in current thermal formulation)

% ==========================================================================
%  SECTION 2.6 – Gas-liquid equilibria & carbonate constants (Eq. 57-63)
% ==========================================================================

T_K    = T + 273.15;   % temperature in Kelvin
T_ref  = p.T_ref;      % 298.15 K
R_gas  = p.R_gas;

% Henry constants (temperature-dependent, Eq. 60)
KH_O2  = p.KH_O2_ref  * exp( p.C_O2  * (1/T_K - 1/T_ref));
KH_CO2 = p.KH_CO2_ref * exp( p.C_CO2 * (1/T_K - 1/T_ref));

% Equilibrium dissolved concentrations (Eq. 57-59)
Xeq_O2  = KH_O2  * p.p_atm * p.y_O2;       % [mol/m³]
COeq_CO2 = KH_CO2 * p.p_atm * p.y_CO2_atm;  % [mol/m³]  atm. equilibrium
COinj_CO2 = KH_CO2 * p.p_atm * p.y_CO2_inj; % [mol/m³]  pure CO2 injection

% Carbonate equilibrium constants (Eq. 61-63)
K1 = p.K1_ref * exp(-p.dH_K1/R_gas * (1/T_K - 1/T_ref));
K2 = p.K2_ref * exp(-p.dH_K2/R_gas * (1/T_K - 1/T_ref));
KW = p.KW_ref * exp(-p.dH_KW/R_gas * (1/T_K - 1/T_ref));

% ==========================================================================
%  SECTION 2.5 – Carbonate speciation (Eq. 40-44)
% ==========================================================================

Delta   = H^2 + H*K1 + K1*K2;                  % (Eq. 40)
CO2_aq  = DIC * H^2 / Delta;                    % (Eq. 41)
HCO3    = DIC * H * K1 / Delta;                 % (Eq. 42)
CO3     = DIC * K1 * K2 / Delta;                % (Eq. 43)
OH      = KW / H;                                % (Eq. 44)

% Ratios used in the proton derivative (Eq. 52)
g_ratio = HCO3 / (DIC + 1e-30);
h_ratio = CO3  / (DIC + 1e-30);

% Partial derivatives of ∆ w.r.t. H (Eq. 55-56)
dDelta_dH  = 2*H + K1;
dg_dH      = K1*(K1*K2 - H^2) / Delta^2;
dh_dH      = -K1*K2 / Delta^2 * dDelta_dH;

% fDIC and fH from electroneutrality (Eq. 50-51)
fDIC = -(g_ratio + 2*h_ratio);
fCat = 1;
fH   = 1 - DIC*(dg_dH + 2*dh_dH) + KW/H^2;

% ==========================================================================
%  SECTION 2.5 – Gas transfer (Eq. 34-39)
% ==========================================================================

% Superficial gas velocities
Ug_CO2 = QCO2 / (p.Asump + 1e-30);
Ug_O2  = Qair  / (p.Asump + 1e-30);

% Volumetric transfer coefficients (power-law, Eq. 34-35)
%   Use smoothmax(x,0) via log-sum-exp to keep symbolic differentiability
eps_reg = 1e-10;
kLa_CO2 = p.alpha_CO2 * (Ug_CO2 + eps_reg)^p.beta_CO2;   % [s⁻¹] in sump
kLa_O2  = p.alpha_O2  * (Ug_O2  + eps_reg)^p.beta_O2;

% Effective coefficients scaled to reactor volume (Eq. 36-37)
kLa_CO2_eff = kLa_CO2 * p.Asump / (V + 1e-30);
kLa_O2_eff  = kLa_O2  * p.Asump / (V + 1e-30);

% Cross-stripping (Eq. 38-39)
k_strip_CO2_eff = p.k_strip_CO2_O2 * kLa_O2_eff;
k_strip_O2_eff  = p.k_strip_O2_CO2 * kLa_CO2_eff;

% Paddlewheel surface exchange (depth-dependent)
Depth_val = (V - p.Vsump) / (p.W * p.L + 1e-30);  % (Eq. 6)

% ==========================================================================
%  SECTION 2.4 – Biological submodel (Eq. 24-31)
% ==========================================================================

% Photosynthetically active radiation (Eq. 2)
PAR = 0.46 * 4.56 * RadG;   % [µmol/(m²·s)]

% Depth-averaged irradiance with Beer-Lambert (Eq. 24)
Iav = PAR / (p.Ka * Depth_val * Xalg + 1e-30) * ...
      (1 - exp(-p.Ka * Depth_val * Xalg));

% Light limitation (Eq. 25)
mu_I = Iav^p.n_hill / (p.Ik^p.n_hill + Iav^p.n_hill + 1e-30);

% Cubic Hermite window function (Eq. 26)
% Generic: given [a, b, c] returns factor in [0,1]
mu_T   = smooth_window(T,   p.T_min,   p.T_opt,   p.T_max);

% Convert H [mol/m³] to pH for biological limitation
pH_sym = -log(H / 1000) / log(10);   % (Eq. 3)
mu_pH  = smooth_window(pH_sym, p.pH_min, p.pH_opt, p.pH_max);

% DO saturation ratio (Eq. 27)
DO_sat = 100 * XO2 / (Xeq_O2 + 1e-30);

% DO inhibition (Eq. 28)
mu_DO  = 1 - (DO_sat / p.DO_max)^p.m_DO;
mu_DO  = fmax(mu_DO, 0);   % clamp to [0,1]

% Gross photosynthetic rate (Eq. 29)
P_rate = p.mu_max * mu_I * mu_T * mu_pH * mu_DO;

% Specific growth rate (Eq. 30)
mu_g   = p.eta_X * P_rate;

% Maintenance / respiration (Eq. 31)
m_resp = p.m_min * (1 + p.k_resp_I * (1 - mu_I)) * ...
         p.Q10^((T - 20) / 10);

% ==========================================================================
%  SECTION 2.3 – Thermal submodel (Eq. 7-23)
% ==========================================================================

% --- Solar irradiance gain (Eq. 9) ---
Q_irrad = p.alpha_rad * p.A * RadG;

% --- Longwave radiative exchange (Eq. 10) ---
%   Tsky from clear-sky emissivity: eps_sky ≈ 0.787 + 0.764*ln(RH/100+ε)·(T_K/273)^(1/7)
%   Simple approximation:
eps_sky = 0.70 + 5.95e-5 * (RH/100 * 4.596 * exp(17.27*Text/(237.3+Text)));
T_sky_K = (eps_sky * (Text+273.15)^4)^(1/4);
Q_rad   = p.sigma * p.eps_w * p.A * (T_sky_K^4 - T_K^4);

% --- Ground conduction (Eq. 11-12) ---
R_pp    = p.x_liner/p.K_liner + p.x_ug/p.K_ug;
Q_cond  = (p.A / R_pp) * (p.T_g - T);

% --- Evaporation (Eq. 13-14) ---
% Vapour concentration at saturation (Antoine approx.) [kg/m³]
e_sat   = 611.2 * exp(17.67 * T / (T + 243.5));        % [Pa]
e_air   = (RH/100) * 611.2 * exp(17.67 * Text / (Text + 243.5));
Cs_vap  = 0.622 * e_sat / ((p.p_atm * 101325) - e_sat) * p.rho_w;
Ca_vap  = 0.622 * e_air / ((p.p_atm * 101325) - e_air) * p.rho_w;
m_dot_e = p.k_m * p.A * fmax(Cs_vap - Ca_vap, 0);      % [kg/s]  (Eq. 13)
lv      = (2.501 - 0.00237 * T) * 1e6;                  % [J/kg]  latent heat
Q_evap  = -lv * m_dot_e;                                % [W]     (Eq. 14)

% Evaporative volume loss
V_dot_e = m_dot_e / p.rho_w;                            % [m³/s]

% --- Convection (Eq. 15) ---
Q_conv = p.h_c * p.A * (Text - T);

% --- Hydraulic enthalpy (Eq. 16-17) ---
Q_dil  =  p.rho_w * p.Cp_w * Qd * p.T_in_medium;       % inflow at T_in
Q_harv = -p.rho_w * p.Cp_w * Qh * T;                   % outflow at T

% --- Paddlewheel mixing heat (Eq. 18) ---
Q_mix  = p.P_mix;

% --- Spiral heat exchanger (Eq. 19-23) ---
Cw_HX  = p.rho_w * p.Cp_w * Qw;                        % (Eq. 21)
eps_HX = 1 - exp(-p.UA_HX / (Cw_HX + 1e-30));          % (Eq. 22)
Q_HX   = Cw_HX * (Tin_hx - T) * eps_HX;                % (Eq. 23)

% Total heat flux
Q_Sigma = Q_irrad + Q_rad + Q_cond + Q_evap + ...
          Q_conv  + Q_dil + Q_harv + Q_mix  + Q_HX;    % (Eq. 8)

% ==========================================================================
%  SECTION 2.5 – ODE RIGHT-HAND SIDE (Eq. 7, 32-47, 49)
% ==========================================================================

% --- Hydraulics (Eq. 32) ---
V_dot  = Qd - Qh - V_dot_e;

% --- Temperature ODE (Eq. 7) ---
T_dot  = Q_Sigma / (p.rho_w * p.Cp_w * V + 1e-30) - ...
         (T / (V + 1e-30)) * V_dot;

% --- Biomass balance (Eq. 33) ---
Xalg_dot = (mu_g - m_resp) * Xalg - (Qd / (V + 1e-30)) * Xalg;

% --- DIC balance (Eq. 45) ---
DIC_dot =  (Qd / (V + 1e-30)) * (p.DIC_in - DIC)                       ...  % dilution
         - P_rate * Xalg / (p.Y_CO2 / p.M_CO2)                          ...  % photosynthesis uptake
         + m_resp * Xalg  / (p.Y_CO2 / p.M_CO2)                         ...  % respiration release
         + kLa_CO2_eff * (COinj_CO2 - CO2_aq)                           ...  % sump CO2 injection (Eq. 45 line 3)
         + p.k_atm_CO2  * (COeq_CO2 - CO2_aq)                           ...  % atmospheric exchange
         + p.k_pw_CO2 * p.W * p.Lpw * Depth_val / (V+1e-30) * (COeq_CO2 - CO2_aq) ...  % paddlewheel
         + k_strip_CO2_eff * (COeq_CO2 - CO2_aq);                             % cross-stripping

% --- Strong cations (Eq. 46) ---
Cat_dot = (Qd / (V + 1e-30)) * (p.Cat_in - Cat);

% --- Dissolved O2 balance (Eq. 47) ---
XO2_dot =  (Qd / (V + 1e-30)) * (Xeq_O2 - XO2)                         ...  % dilution at sat.
          + P_rate * Xalg  / (p.Y_O2 / p.M_O2)                          ...  % photosynthetic production
          - m_resp * Xalg  / (p.Y_O2 / p.M_O2)                          ...  % respiration consumption
          + kLa_O2_eff  * (Xeq_O2 - XO2)                                ...  % sump aeration
          + p.k_atm_O2  * (Xeq_O2 - XO2)                                ...  % atmospheric
          + p.k_pw_O2 * p.W * p.Lpw * Depth_val / (V+1e-30) * (Xeq_O2 - XO2) ...  % paddlewheel
          - k_strip_O2_eff * XO2;                                              % cross-stripping loss

% --- Proton dynamics from electroneutrality (Eq. 49) ---
H_dot  = -(fDIC * DIC_dot + fCat * Cat_dot) / (fH + 1e-30);

% --- Assemble state derivative ---
xdot = vertcat(Xalg_dot, XO2_dot, DIC_dot, Cat_dot, H_dot, T_dot, V_dot);

% ==========================================================================
%  OUTPUT MAP  y = [pH, DO%, Xalg_gL, Depth]  (Eq. 3-6)
% ==========================================================================

pH_out    = -log(H / 1000) / log(10);                   % (Eq. 3)
DO_out    = 100 * XO2 / (Xeq_O2 + 1e-30);              % (Eq. 4)
Xalg_gL  = Xalg / 1000;                                 % (Eq. 5)
Depth_out = (V - p.Vsump) / (p.W * p.L + 1e-30);       % (Eq. 6)

y = vertcat(pH_out, DO_out, Xalg_gL, Depth_out);

% ==========================================================================
%  CASADI FUNCTION OBJECTS
% ==========================================================================

ode_fun = Function('ode_fun', {x, u, d}, {xdot}, ...
                   {'x','u','d'}, {'xdot'});

out_fun = Function('out_fun', {x, u, d}, {y}, ...
                   {'x','u','d'}, {'y'});

end % build_raceway_ode

%% =========================================================================
%  LOCAL HELPER: smooth_window
%  Cubic Hermite window function (Eq. 26) for temperature/pH limitation.
%  Returns µ ∈ [0,1] given variable val and window [a, b, c].
% =========================================================================

function mu = smooth_window(val, a, b, c)
% Cubic Hermite: µ = 3r² - 2r³, with appropriate r in each sub-interval

import casadi.*

% r for left branch  (a < val <= b):  r = (val-a)/(b-a)
r_left  = (val - a) / (b - a + 1e-30);
mu_left = 3*r_left^2 - 2*r_left^3;

% r for right branch (b < val < c):   r = (c-val)/(c-b)
r_right  = (c - val) / (c - b + 1e-30);
mu_right = 3*r_right^2 - 2*r_right^3;

% Select branch smoothly via if_else; outside [a,c] → 0
mu_in_left  = if_else(val > a,  mu_left,  0);
mu_at_left  = if_else(val <= b, mu_in_left, 0);
mu_in_right = if_else(val > b,  mu_right, 0);
mu_at_right = if_else(val < c,  mu_in_right, 0);

mu = mu_at_left + mu_at_right;

% Clamp to [0, 1] for robustness
mu = fmin(fmax(mu, 0), 1);

end % smooth_window
