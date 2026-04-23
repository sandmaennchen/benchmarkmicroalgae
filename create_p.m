function p = create_p()
p = struct();
% Geometry
if ~isfield(p,'L'), p.L = 80; end                              % [m]         channel length
if ~isfield(p,'W'), p.W = 1; end                               % [m]         channel width
if ~isfield(p,'A'), p.A = p.L * p.W; end                      % [m²]        illuminated surface
if ~isfield(p,'rsump'), p.rsump = 0.4; end                     % [m]         sump radius
if ~isfield(p,'hsump'), p.hsump = 1.12; end                    % [m]         sump depth
if ~isfield(p,'Asump'), p.Asump = pi * p.rsump^2; end          % [m²]        sump cross-section
if ~isfield(p,'Vsump'), p.Vsump = p.Asump * p.hsump; end       % [m³]        sump volume

% Water and thermal
if ~isfield(p,'rho_w'), p.rho_w = 997; end                     % [kg/m³]     water density
if ~isfield(p,'Cp_w'), p.Cp_w = 4182; end                      % [J/(kg·K)]  water heat capacity
if ~isfield(p,'alpha_rad'), p.alpha_rad = 0.97; end             % [-]         shortwave absorptance
if ~isfield(p,'sigma'), p.sigma = 5.6697e-8; end               % [W/(m²·K⁴)] Stefan-Boltzmann constant
if ~isfield(p,'eps_w'), p.eps_w = 0.9875; end                  % [-]         water emissivity
if ~isfield(p,'x_liner'), p.x_liner = 0.002; end               % [m]         liner thickness
if ~isfield(p,'K_liner'), p.K_liner = 0.33; end                % [W/(m·K)]   liner conductivity
if ~isfield(p,'x_ug'), p.x_ug = 2.3; end                       % v [m]         subgrade thickness
if ~isfield(p,'K_ug'), p.K_ug = 4.45; end                      % v [W/(m·K)]   subgrade conductivity
if ~isfield(p,'T_g'), p.T_g = 17.4; end                        % [°C]        ground reference temperature
if ~isfield(p,'k_m'), p.k_m = 0.01; end                        % [m/s]       evaporation mass-transfer coeff.
if ~isfield(p,'h_c'), p.h_c = 10; end                          % [W/(m²·K)]  convective heat-transfer coeff.
if ~isfield(p,'P_mix'), p.P_mix = 3389; end                    % v [W]         paddlewheel mechanical power
if ~isfield(p,'T_in_medium'), p.T_in_medium = 20; end          % [°C]        dilution medium temperature
if ~isfield(p,'Lpw'), p.Lpw = 2; end                           % [m]         paddlewheel zone effective length

% Biology
if ~isfield(p,'mu_max'), p.mu_max = 2.23e-5; end             % v [s⁻¹]       max. specific photosynthesis (0.8 d⁻¹)
if ~isfield(p,'Ka'), p.Ka = 0.0824; end                           % v [m²/g]      light extinction coefficient
if ~isfield(p,'Ik'), p.Ik = 122.27; end                           % v [µmol/(m²·s)] half-sat. irradiance
if ~isfield(p,'n_hill'), p.n_hill = 2.86; end                      % v [-]         Hill exponent
if ~isfield(p,'T_min'), p.T_min = -0.9; end                    % v [°C]        temperature window minimum
if ~isfield(p,'T_opt'), p.T_opt = 30.72; end                   % v [°C]        optimal temperature
if ~isfield(p,'T_max'), p.T_max = 40.81; end                   % v [°C]        temperature window maximum
if ~isfield(p,'pH_min'), p.pH_min = 5.33; end                  % v [-]         pH window minimum
if ~isfield(p,'pH_opt'), p.pH_opt = 7.8; end                   % v [-]         optimal pH
if ~isfield(p,'pH_max'), p.pH_max = 10.39; end                 % v [-]         pH window maximum
if ~isfield(p,'DO_max'), p.DO_max = 304.15; end                % v [%]         DO inhibition threshold
if ~isfield(p,'m_DO'), p.m_DO = 1.54; end                         % v [-]         inhibition exponent
if ~isfield(p,'eta_X'), p.eta_X = 1.0; end                     % v [-]         gross-to-growth allocation
if ~isfield(p,'Y_CO2'), p.Y_CO2 = 44/30; end                 % [gCO2/gXalg] CO2 yield coefficient
if ~isfield(p,'Y_O2'), p.Y_O2 = 32/30; end                   % [gO2/gXalg]  O2 yield coefficient
if ~isfield(p,'M_CO2'), p.M_CO2 = 44; end                      % [g/mol]     molar mass of CO2
if ~isfield(p,'M_O2'), p.M_O2 = 32; end                       % [g/mol]     molar mass of O2
if ~isfield(p,'m_min'), p.m_min = 6.6e-7; end              % v [s⁻¹]       basal respiration rate at 20°C
if ~isfield(p,'k_resp_I'), p.k_resp_I = 2.6e-7; end               % v [-]         low-light maintenance gain
if ~isfield(p,'Q10'), p.Q10 = 2.17; end                         % v [-]         temperature sensitivity (Q10)

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
if ~isfield(p,'dH_K1'), p.dH_K1 = -7.7e3; end                   % [J/mol]     enthalpy of K1
if ~isfield(p,'dH_K2'), p.dH_K2 = -14.9e3; end                 % [J/mol]     enthalpy of K2
if ~isfield(p,'dH_KW'), p.dH_KW = 55.8e3; end                  % [J/mol]     enthalpy of KW
if ~isfield(p,'alpha_CO2'), p.alpha_CO2 = 3.75; end            % v [-]         kLa scale factor for CO2
if ~isfield(p,'beta_CO2'), p.beta_CO2 = 1.12; end               % v [-]         kLa power-law exponent for CO2
if ~isfield(p,'alpha_O2'), p.alpha_O2 = 1.36; end              % v [-]         kLa scale factor for O2
if ~isfield(p,'beta_O2'), p.beta_O2 = 0.95; end                 % v [-]         kLa power-law exponent for O2
if ~isfield(p,'k_strip_CO2_O2'), p.k_strip_CO2_O2 = 0.25; end  % v [-]         CO2 cross-stripping factor
if ~isfield(p,'k_strip_O2_CO2'), p.k_strip_O2_CO2 = 0.14; end  % v [-]         O2 cross-stripping factor
if ~isfield(p,'k_atm_CO2'), p.k_atm_CO2 = 1.35e-4; end         % v [s⁻¹]       atmospheric CO2 exchange coeff.
if ~isfield(p,'k_atm_O2'), p.k_atm_O2 = 2.24e-4; end           % v [s⁻¹]       atmospheric O2 exchange coeff.
if ~isfield(p,'k_pw_CO2'), p.k_pw_CO2 = 5.38e-4; end           % v [s⁻¹]       paddlewheel CO2 exchange coeff.
if ~isfield(p,'k_pw_O2'), p.k_pw_O2 = 3.48e-4; end             % v [s⁻¹]       paddlewheel O2 exchange coeff.
if ~isfield(p,'DIC_in'), p.DIC_in = 0.017*1e3; end                     % vt [mol/m³]    inlet dissolved inorganic carbon
if ~isfield(p,'Cat_in'), p.Cat_in = 0.011*1e3; end                     % vt [mol/m³]    inlet strong cations
% additional
if ~isfield(p,'reg_eps'), p.reg_eps = 1e-12; end                     % [-]    regularization
if ~isfield(p,'R_w'), p.R_w = 9.375e-5; end                     % [m^2 * K * W^{-1}]    HX wall conduction resistance
if ~isfield(p,'R_o'), p.R_o = 0.0013; end                     % [m^2 * K * W^{-1}]    outer convective resistance
if ~isfield(p,'R_i'), p.R_i = 0.0013; end                     % ???? [m^2 * K * W^{-1}]    inner convective resistance
if ~isfield(p,'A_o'), p.A_o = 2.3876; end                     % [m^2] HX external area

% TODO: in table but not used yet
if ~isfield(p,'c_evap'), p.c_evap = 0.89; end                     % v [-] evaporation empirical constant
if ~isfield(p,'c_wind'), p.c_wind = 0.93; end                     % v [-] effective wind constant

p.U = 1/(p.R_i + p.R_w + p.R_o);
p.UA = p.U * p.A_o;
p.UA_HX = p.UA; % not v
end