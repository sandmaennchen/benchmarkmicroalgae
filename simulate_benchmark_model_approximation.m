function [results] = simulate_benchmark_model_approximation(Data, ctrl)
% Approximate benchmark simulator with CasADi-integrated continuous dynamics.
% Uses build_ode(p) inside the closed-loop simulation at each sample.

[t, Env] = build_time_and_env(Data);
N = numel(t);

dt = 60;

refs = struct('pH', 7.5, 'DO', 150.0, 'T', 25.0);

p = create_p();

if ~isempty(Env.Temp_ext)
    p.T_in_medium = Env.Temp_ext(1);
end

[ode_fun, out_fun, growth_terms_fun] = build_ode(p);
F_step = build_integrator(ode_fun, dt);

x0 = get_initial_state(p, Env);
xk = x0;

pH = zeros(N,1);
DO = zeros(N,1);
Depth = zeros(N,1);
X_gL = zeros(N,1);

Xalg = zeros(N,1);
XO2 = zeros(N,1);
DIC = zeros(N,1);
Cat = zeros(N,1);
H = zeros(N,1);
T = zeros(N,1);
V = zeros(N,1);

mu_I = zeros(N, 1);
mu_T = zeros(N, 1);
mu_pH = zeros(N, 1);
mu_DO = zeros(N, 1);
mu = zeros(N, 1);
m = zeros(N, 1);
P = zeros(N, 1);

QCO2_cmd = zeros(N,1);
Qair_cmd = zeros(N,1);
QCO2_del = zeros(N,1);
Qair_del = zeros(N,1);
Qd_bin   = zeros(N,1);
Qh_bin   = zeros(N,1);
Qd       = zeros(N,1);
Qh       = zeros(N,1);
Qhx_cmd  = zeros(N,1);
Tin_hx   = zeros(N,1);

st_CtrlSignals = struct();
state_pH_CO2 = struct();
state_DO_air = struct();
state_HD = struct();
state_Temp_HX = struct();

% QCO2_max = 20/1000/60;
% Qair_max = 500/1000/60;
Qd_rate  = 20/1000/60;
Qh_rate  = 20/1000/60;
Qw_min = 0;
Qw_max = 1e-2;
Tin_min = 0;
Tin_max = 50;

for k = 1:N
    tk = t(k);
    secday = rem(tk, 86400);
    Timeline = struct();
    Timeline.dt = dt;
    Timeline.index = k;
    Timeline.time = tk;
    Timeline.time_secday = secday;
    Timeline.hour = floor(secday/3600);
    Timeline.min = floor(rem(secday,3600)/60);

    env = struct();
    env.RadGlobal = Env.RadG(k);
    env.RadPAR = Env.RadPAR(k);
    env.Temp_ext = Env.Temp_ext(k);
    env.RH = Env.RH(k);
    env.Wind = Env.Wind(k);

    dk = [Env.RadG(k); Env.RH(k); Env.Temp_ext(k); Env.Wind(k)];
    yk = full(out_fun(xk, zeros(6,1), dk));
    pH(k) = yk(1);
    DO(k) = yk(2);
    X_gL(k) = yk(3);
    Depth(k) = yk(4);

    T(k) = xk(6);
    
    obs = struct();
    obs.pH = pH(k);
    obs.DO = DO(k);
    obs.Depth = Depth(k);
    obs.Xalg_gL = X_gL(k);
    obs.T = T(k);
    % obs.T = Env.Temp_ext(k);
    % obs

    Xalg(k) = xk(1);
    XO2(k) = xk(2);
    DIC(k) = xk(3);
    Cat(k) = xk(4);
    H(k) = xk(5);
    T(k) = xk(6);
    V(k) = xk(7);


    future = project_future(t, Env, k);

    [st_CtrlSignals, state_pH_CO2] = ctrl.fn_pH_CO2(Timeline, obs, refs, env, future, st_CtrlSignals, state_pH_CO2); % Qco2

    [st_CtrlSignals, state_DO_air] = ctrl.fn_DO_air(Timeline, obs, refs, env, future, st_CtrlSignals, state_DO_air); % Qair

    [st_CtrlSignals, state_HD] = ctrl.fn_HD(Timeline, obs, refs, env, future, st_CtrlSignals, state_HD); % Qd_bin, Qh_bin

    [st_CtrlSignals, state_Temp_HX] = ctrl.fn_Temp_HX(Timeline, obs, refs, env, future, st_CtrlSignals, state_Temp_HX); % Qhx, Tin_hx

    QCO2_cmd(k) = st_CtrlSignals.Qco2;
    Qair_cmd(k) = st_CtrlSignals.Qair;
    Qd_bin(k)   = st_CtrlSignals.Qd_bin;
    Qh_bin(k)   = st_CtrlSignals.Qh_bin;
    Qhx_cmd(k)  = st_CtrlSignals.Qhx;
    Tin_hx(k)   = st_CtrlSignals.Tin_hx;

    QCO2_del(k) = QCO2_cmd(k);
    Qair_del(k) = Qair_cmd(k);
    Qd(k) = Qd_bin(k) * Qd_rate;
    Qh(k) = Qh_bin(k) * Qh_rate;

    uk = [QCO2_del(k); Qair_del(k); Qd(k); Qh(k); Qhx_cmd(k); Tin_hx(k)];

    gt = full(growth_terms_fun(xk, uk, dk));
    mu_I(k) = gt(1);
    mu_T(k) = gt(2);
    mu_pH(k) = gt(3);
    mu_DO(k) = gt(4);
    mu(k) = gt(5);
    m(k) = gt(6);
    P(k) = gt(7);
    step_out = F_step('x0', xk, 'p', [uk; dk]);
    xk = full(step_out.xf);
end

cum_air_L = cumsum(Qair_del * dt * 1000);
cum_CO2_L = cumsum(QCO2_del * dt * 1000);

area_m2 = p.A;
days = (t(end) - t(1)) / 86400;

V_L = Depth * area_m2 * 1000;

mass_in_pond_start = X_gL(1) * V_L(1);
mass_in_pond_end   = X_gL(end) * V_L(end);

gain_in_pond_g = mass_in_pond_end - mass_in_pond_start;
cum_harv_g = cumsum(Qh * dt * 1000 .* X_gL);
harv_total_g = cum_harv_g(end);
prod_total_g = gain_in_pond_g + harv_total_g;
prod_areal_gm2_day = prod_total_g / (days * area_m2);
harv_frac = 100 * (harv_total_g / prod_total_g);
acumm_rel = 100 * (X_gL(end) - X_gL(1)) / X_gL(1);

harv_prod_areal_gm2_day = harv_total_g / max(days * area_m2, eps);

[CO2, HCO3, CO3] = i_carbonate_species(DIC, H, T, p);

results = struct();
results.t = t;
results.refs = refs;

results.pH = pH;
results.DO = DO;
results.T = T;
results.X_gL = X_gL;
results.Depth = Depth;

results.QCO2_cmd = QCO2_cmd;
results.Qair_cmd = Qair_cmd;
results.QCO2_del = QCO2_del;
results.Qair_del = Qair_del;
results.Qd = Qd;
results.Qh = Qh;

results.cum_air_L = cum_air_L;
results.cum_CO2_L = cum_CO2_L;
results.cum_harv_g = cum_harv_g;
results.total_air_L = cum_air_L(end);
results.total_CO2_L = cum_CO2_L(end);

results.gain_g = gain_in_pond_g;
results.acumm_rel = acumm_rel;
results.prod_g = prod_total_g;
results.prod_areal_gm2_day = prod_areal_gm2_day;
results.harv_total_g = harv_total_g;
results.harv_frac = harv_frac;
results.harv_prod_areal_gm2_day = harv_prod_areal_gm2_day;

results.mu_I = mu_I;
results.mu_T = mu_T;
results.mu_pH = mu_pH;
results.mu_DO = mu_DO;
results.P = P;
results.mu = mu;
results.m = m;

results.DIC = DIC;
results.Cat = Cat;
results.HCO3 = HCO3;
results.CO3 = CO3;
results.CO2 = CO2;

results.HX = struct();
results.HX.Qw_m3s = min(max(Qhx_cmd, Qw_min), Qw_max);
results.HX.Tin_C = min(max(Tin_hx, Tin_min), Tin_max);
results.HX.Tout_C = T;
results.HX.UA_WK = (p.UA_HX/1000) * ones(N,1);
results.HX.Q_W = p.rho_w * p.Cp_w .* results.HX.Qw_m3s .* (results.HX.Tin_C - T);
results.HX.limits = struct('Qw_min',Qw_min,'Qw_max',Qw_max,'Tin_min',Tin_min,'Tin_max',Tin_max);

results.Env = Env;

results.J = struct();
results.J.pH = mean((pH - refs.pH).^2);
results.J.DO = mean((DO - refs.DO).^2);
results.J.T  = mean((T - refs.T).^2);


end

function [CO2, HCO3, CO3] = i_carbonate_species(DIC, H, T, p)
K1 = p.K1_ref * exp(-p.dH_K1/p.R_gas * (1./(T+273.15) - 1/p.T_ref));
K2 = p.K2_ref * exp(-p.dH_K2/p.R_gas * (1./(T+273.15) - 1/p.T_ref));

Delta = H.^2 + H.*K1 + K1.*K2;
CO2 = DIC .* H.^2 ./ max(Delta, 1e-16);
HCO3 = DIC .* H .* K1 ./ max(Delta, 1e-16);
CO3 = DIC .* K1 .* K2 ./ max(Delta, 1e-16);

end
