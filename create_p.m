function p = create_p()
p = struct();
% Geometry
p.L = 80;                              % [m]         channel length
p.W = 1;                               % [m]         channel width
p.A = p.L * p.W;                      % [m²]        illuminated surface
p.rsump = 0.4;                     % [m]         sump radius
p.hsump = 1.12;                    % [m]         sump depth
p.Asump = pi * p.rsump^2;          % [m²]        sump cross-section
p.Vsump = p.Asump * p.hsump;       % [m³]        sump volume

% Water and thermal
p.rho_w = 997;                     % [kg/m³]     water density
p.Cp_w = 4182;                      % [J/(kg·K)]  water heat capacity
p.alpha_rad = 0.97;             % [-]         shortwave absorptance
p.sigma = 5.6697e-8;               % [W/(m²·K⁴)] Stefan-Boltzmann constant
p.eps_w = 0.9875;                  % [-]         water emissivity
p.x_liner = 0.002;               % [m]         liner thickness
p.K_liner = 0.33;                % [W/(m·K)]   liner conductivity
p.x_ug = 2.3;                       % v [m]         subgrade thickness
p.K_ug = 4.45;                      % v [W/(m·K)]   subgrade conductivity
p.T_g = 17.4;                        % [°C]        ground reference temperature
p.k_m = 0.01;                        % [m/s]       evaporation mass-transfer coeff.
p.h_c = 10;                          % [W/(m²·K)]  convective heat-transfer coeff.
p.P_mix = 3389;                    % v [W]         paddlewheel mechanical power
p.T_in_medium = 20;          % [°C]        dilution medium temperature
p.Lpw = 2;                           % [m]         paddlewheel zone effective length

% Biology
p.mu_max = 2.23e-5;             % v [s⁻¹]       max. specific photosynthesis (0.8 d⁻¹)
p.Ka = 0.0824;                           % v [m²/g]      light extinction coefficient
p.Ik = 122.27;                           % v [µmol/(m²·s)] half-sat. irradiance
p.n_hill = 2.86;                      % v [-]         Hill exponent
p.T_min = -0.9;                    % v [°C]        temperature window minimum
p.T_opt = 30.72;                   % v [°C]        optimal temperature
p.T_max = 40.81;                   % v [°C]        temperature window maximum
p.pH_min = 5.33;                  % v [-]         pH window minimum
p.pH_opt = 7.8;                   % v [-]         optimal pH
p.pH_max = 10.39;                 % v [-]         pH window maximum
p.DO_max = 304.15;                % v [%]         DO inhibition threshold
p.m_DO = 1.54;                         % v [-]         inhibition exponent
p.eta_X = 1.0;                     % v [-]         gross-to-growth allocation
p.Y_CO2 = 44/30;                 % [gCO2/gXalg] CO2 yield coefficient
p.Y_O2 = 32/30;                   % [gO2/gXalg]  O2 yield coefficient
p.M_CO2 = 44;                      % [g/mol]     molar mass of CO2
p.M_O2 = 32;                       % [g/mol]     molar mass of O2
p.m_min = 6.6e-7;              % v [s⁻¹]       basal respiration rate at 20°C
p.k_resp_I = 2.6e-7;               % v [-]         low-light maintenance gain
p.Q10 = 2.17;                         % v [-]         temperature sensitivity (Q10)

% Chemistry and transfer
p.T_ref = 298.15;                  % [K]         reference temperature
p.R_gas = 8.314;                   % [J/(mol·K)] ideal gas constant
p.p_atm = 1.0;                     % [atm]       atmospheric pressure
p.y_O2 = 0.2095;                      % [-]         O2 molar fraction in air
p.y_CO2_atm = 0.00039;         % [-]         CO2 molar fraction in air (~380 ppm)
p.y_CO2_inj = 1.0;            % [-]         CO2 molar fraction in injection gas (pure CO2)
p.KH_O2_ref = 1.3e-3 * 1e3;   % [mol/(m³·atm)] Henry const. for O2 at T_ref
p.KH_CO2_ref = 3.4e-2 * 1e3; % [mol/(m³·atm)] Henry const. for CO2 at T_ref
p.C_O2 = 1700;                      % [K]         Van't Hoff constant for O2
p.C_CO2 = 2400;                   % [K]         Van't Hoff constant for CO2
p.K1_ref = 10^(-6.3) * 1e3;         % [mol/m³]    first carbonate dissociation const.
p.K2_ref = 10^(-10.33) * 1e3;        % [mol/m³]    second carbonate dissociation const.
p.KW_ref = 1e-14 * 1e6;          % [mol²/m⁶]   water autoprotolysis const.
p.dH_K1 = -7.7e3;                   % [J/mol]     enthalpy of K1
p.dH_K2 = -14.9e3;                 % [J/mol]     enthalpy of K2
p.dH_KW = 55.8e3;                  % [J/mol]     enthalpy of KW
p.alpha_CO2 = 3.75;            % v [-]         kLa scale factor for CO2
p.beta_CO2 = 1.12;               % v [-]         kLa power-law exponent for CO2
p.alpha_O2 = 1.36;              % v [-]         kLa scale factor for O2
p.beta_O2 = 0.95;                 % v [-]         kLa power-law exponent for O2
p.k_strip_CO2_O2 = 0.25;  % v [-]         CO2 cross-stripping factor
p.k_strip_O2_CO2 = 0.14;  % v [-]         O2 cross-stripping factor
p.k_atm_CO2 = 1.35e-4;         % v [s⁻¹]       atmospheric CO2 exchange coeff.
p.k_atm_O2 = 2.24e-4;           % v [s⁻¹]       atmospheric O2 exchange coeff.
p.k_pw_CO2 = 5.38e-4;           % v [s⁻¹]       paddlewheel CO2 exchange coeff.
p.k_pw_O2 = 3.48e-4;             % v [s⁻¹]       paddlewheel O2 exchange coeff.
p.DIC_in = 0.017*1e3;                     % vt [mol/m³]    inlet dissolved inorganic carbon
p.Cat_in = 0.011*1e3;                     % vt [mol/m³]    inlet strong cations
% additional
p.reg_eps = 1e-12;                     % [-]    regularization
p.R_w = 9.375e-5;                     % [m^2 * K * W^{-1}]    HX wall conduction resistance
p.R_o = 0.0013;                     % [m^2 * K * W^{-1}]    outer convective resistance
p.R_i = 0.0013;                     % ???? [m^2 * K * W^{-1}]    inner convective resistance
p.A_o = 2.3876;                     % [m^2] HX external area

% TODO: in table but not used yet
p.c_evap = 0.89;                     % v [-] evaporation empirical constant
p.c_wind = 0.93;                     % v [-] effective wind constant

p.U = 1/(p.R_i + p.R_w + p.R_o);
p.UA = p.U * p.A_o;
p.UA_HX = p.UA; % not v
end