function [results] = simulate_benchmark_model_approximation(Data, ctrl)
% Approximate benchmark simulator with CasADi-integrated continuous dynamics.
% Uses build_ode(p) inside the closed-loop simulation at each sample.

if nargin < 2
    error('simulate_benchmark_model_approximation requires Data and ctrl inputs.');
end

if exist('casadi.MX', 'class') ~= 8
    error(['CasADi is not available in MATLAB path. ' ...
           'Install CasADi and ensure import works before running this simulator.']);
end

[t, Env] = i_build_time_and_env(Data);
N = numel(t);
if N < 1
    error('Data does not contain valid samples.');
end

if N > 1
    dt = max(1, t(2) - t(1));
else
    dt = 60;
end

refs = struct('pH', 8.0, 'DO', 150.0, 'T', 30.0);

p = i_default_params();
if ~isempty(Env.Temp_ext)
    p.T_in_medium = Env.Temp_ext(1);
end

[ode_fun, out_fun] = build_ode(p);
F_step = i_build_integrator(ode_fun, dt);

xk = i_initial_state(p, Env);

pH = zeros(N,1);
DO = zeros(N,1);
Depth = zeros(N,1);
X_gL = zeros(N,1);
T = zeros(N,1);

Xalg = zeros(N,1);
XO2 = zeros(N,1);
DIC = zeros(N,1);
Cat = zeros(N,1);
H = zeros(N,1);
V = zeros(N,1);

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

st_pH = struct();
st_DO = struct();
st_HD = struct();
st_TX = struct();

QCO2_max = 20/1000/60;
Qair_max = 500/1000/60;
Qd_rate  = 20/1000/60;
Qh_rate  = 20/1000/60;
Qw_min = 0;
Qw_max = 1e-2;
Tin_min = 0;
Tin_max = 50;

for k = 1:N
    tk = t(k);
    secday = rem(tk, 86400);

    dk = [Env.RadG(k); Env.RH(k); Env.Temp_ext(k); Env.Wind(k)];
    yk = full(out_fun(xk, zeros(6,1), dk));

    pH(k) = yk(1);
    DO(k) = yk(2);
    X_gL(k) = yk(3);
    Depth(k) = yk(4);
    T(k) = xk(6);

    Xalg(k) = xk(1);
    XO2(k) = xk(2);
    DIC(k) = xk(3);
    Cat(k) = xk(4);
    H(k) = xk(5);
    V(k) = xk(7);

    Timeline = struct();
    Timeline.dt = dt;
    Timeline.index = k;
    Timeline.time = tk;
    Timeline.time_secday = secday;
    Timeline.hour = floor(secday/3600);
    Timeline.min = floor(rem(secday,3600)/60);

    obs = struct();
    obs.pH = pH(k);
    obs.DO = DO(k);
    obs.Depth = Depth(k);
    obs.Xalg_gL = X_gL(k);
    obs.T = T(k);

    env = struct();
    env.RadGlobal = Env.RadG(k);
    env.RadPAR = Env.RadPAR(k);
    env.Temp_ext = Env.Temp_ext(k);
    env.RH = Env.RH(k);
    env.Wind = Env.Wind(k);

    future = i_build_future(t, Env, k);

    st_cmd = struct('Qco2',0,'Qair',0,'Qd_bin',0,'Qh_bin',0,'Qhx',0,'Tin_hx',obs.T);

    if isfield(ctrl, 'fn_pH_CO2') && ~isempty(ctrl.fn_pH_CO2)
        [st_tmp, st_pH] = ctrl.fn_pH_CO2(Timeline, obs, refs, env, future, st_cmd, st_pH);
        st_cmd.Qco2 = i_getfield_or(st_tmp, 'Qco2', st_cmd.Qco2);
    end
    if isfield(ctrl, 'fn_DO_air') && ~isempty(ctrl.fn_DO_air)
        [st_tmp, st_DO] = ctrl.fn_DO_air(Timeline, obs, refs, env, future, st_cmd, st_DO);
        st_cmd.Qair = i_getfield_or(st_tmp, 'Qair', st_cmd.Qair);
    end
    if isfield(ctrl, 'fn_HD') && ~isempty(ctrl.fn_HD)
        [st_tmp, st_HD] = ctrl.fn_HD(Timeline, obs, refs, env, future, st_cmd, st_HD);
        st_cmd.Qd_bin = i_getfield_or(st_tmp, 'Qd_bin', st_cmd.Qd_bin);
        st_cmd.Qh_bin = i_getfield_or(st_tmp, 'Qh_bin', st_cmd.Qh_bin);
    end
    if isfield(ctrl, 'fn_Temp_HX') && ~isempty(ctrl.fn_Temp_HX)
        [st_tmp, st_TX] = ctrl.fn_Temp_HX(Timeline, obs, refs, env, future, st_cmd, st_TX);
        st_cmd.Qhx = i_getfield_or(st_tmp, 'Qhx', st_cmd.Qhx);
        st_cmd.Tin_hx = i_getfield_or(st_tmp, 'Tin_hx', st_cmd.Tin_hx);
    end

    QCO2_cmd(k) = i_getfield_or(st_cmd, 'Qco2', 0);
    Qair_cmd(k) = i_getfield_or(st_cmd, 'Qair', 0);
    Qd_bin(k)   = i_getfield_or(st_cmd, 'Qd_bin', 0);
    Qh_bin(k)   = i_getfield_or(st_cmd, 'Qh_bin', 0);
    Qhx_cmd(k)  = i_getfield_or(st_cmd, 'Qhx', 0);
    Tin_hx(k)   = i_getfield_or(st_cmd, 'Tin_hx', obs.T);

    QCO2_del(k) = min(max(QCO2_cmd(k), 0), QCO2_max);
    Qair_del(k) = min(max(Qair_cmd(k), 0), Qair_max);

    if isfield(st_cmd, 'Qd') && ~isempty(st_cmd.Qd)
        Qd(k) = max(0, st_cmd.Qd);
    else
        Qd(k) = Qd_rate * (Qd_bin(k) > 0.5);
    end

    if isfield(st_cmd, 'Qh') && ~isempty(st_cmd.Qh)
        Qh(k) = max(0, st_cmd.Qh);
    else
        Qh(k) = Qh_rate * (Qh_bin(k) > 0.5);
    end

    Qw_k = min(max(Qhx_cmd(k), Qw_min), Qw_max);
    Tin_k = min(max(Tin_hx(k), Tin_min), Tin_max);

    uk = [QCO2_del(k); Qair_del(k); Qd(k); Qh(k); Qw_k; Tin_k];

    try
        step_out = F_step('x0', xk, 'p', [uk; dk]);
        xk = full(step_out.xf);
    catch ME
        error('CasADi integration failed at k=%d: %s', k, ME.message);
    end

    xk = i_state_clip(xk, p);
end

cum_air_L = cumsum(Qair_del * dt * 1000);
cum_CO2_L = cumsum(QCO2_del * dt * 1000);

area_m2 = p.A;
days = max(dt, t(end) - t(1)) / 86400;

V_L = Depth * area_m2 * 1000;
gain_g = (X_gL(end) * V_L(end)) - (X_gL(1) * V_L(1));
prod_g = max(gain_g, 0);
cum_harv_g = cumsum(Qh * dt * 1000 .* X_gL);
harv_total_g = cum_harv_g(end);

acumm_rel = 100 * (X_gL(end) - X_gL(1)) / max(X_gL(1), eps);
prod_areal_gm2_day = prod_g / max(days * area_m2, eps);
harv_prod_areal_gm2_day = harv_total_g / max(days * area_m2, eps);
harv_frac = 100 * harv_total_g / max(prod_g, eps);

[mu_I, mu_T, mu_pH, mu_DO, mu, m, P] = i_growth_terms(Xalg, XO2, DIC, H, T, Depth, Env, p);
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

results.gain_g = gain_g;
results.acumm_rel = acumm_rel;
results.prod_g = prod_g;
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

function F_step = i_build_integrator(ode_fun, dt)
import casadi.*

x = MX.sym('x', 7);
up = MX.sym('up', 10);
u = up(1:6);
d = up(7:10);
xdot = ode_fun(x, u, d);

dae = struct('x', x, 'p', up, 'ode', xdot);
opts = struct('tf', dt, 'abstol', 1e-8, 'reltol', 1e-6);
F_step = integrator('F_step', 'cvodes', dae, opts);

end

function x0 = i_initial_state(p, Env)
Xalg0 = 0.5 * 1e3; % tv: ref results 0.5; TABLE: 0.32 [g/l]
XO20 = p.KH_O2_ref * p.p_atm * p.y_O2;
DIC0 = 0.014; % v
Cat0 = 1.672669964353326e+01; % p.Cat_in;
H0 = 1.582e-05; % v
T0 = Env.Temp_ext(1);
Depth0 = 0.15;
V0 = max(p.Vsump + 1e-2, p.A * Depth0 + p.Vsump);
x0 = [Xalg0; XO20; DIC0; Cat0; H0; T0; V0];

end

function x = i_state_clip(x, p)
x(1) = max(x(1), 1e-8);
x(2) = max(x(2), 1e-10);
x(3) = max(x(3), 1e-10);
x(4) = max(x(4), 1e-10);
x(5) = min(max(x(5), 1e-12), 1);
x(6) = min(max(x(6), -5), 60);
x(7) = max(x(7), p.Vsump + 1e-3);

end

function p = i_default_params()
p = struct();
p.L = 40;
p.W = 1;
p.A = p.L * p.W;
p.rsump = 0.4;
p.hsump = 1.12;
p.Asump = pi * p.rsump^2;
p.Vsump = p.Asump * p.hsump;

p.rho_w = 1000;
p.Cp_w = 4182;
p.alpha_rad = 0.7;
p.sigma = 5.67e-8;
p.eps_w = 0.97;
p.x_liner = 0.003;
p.K_liner = 0.35;
p.x_ug = 0.20;
p.K_ug = 1.0;
p.T_g = 20;
p.k_m = 0.01;
p.h_c = 10;
p.P_mix = 200;
p.UA_HX = 500;
p.T_in_medium = 20;
p.Lpw = 2;

p.mu_max = 0.8/86400;
p.Ka = 0.2;
p.Ik = 200;
p.n_hill = 1;
p.T_min = 5;
p.T_opt = 30;
p.T_max = 40;
p.pH_min = 6;
p.pH_opt = 8;
p.pH_max = 10;
p.DO_max = 383.21;
p.m_DO = 2;
p.eta_X = 0.6;
p.Y_CO2 = 1.92;
p.Y_O2 = 1.41;
p.M_CO2 = 44;
p.M_O2 = 32;
p.m_min = 0.02/86400;
p.k_resp_I = 0.5;
p.Q10 = 2.0;

p.T_ref = 298.15;
p.R_gas = 8.314;
p.p_atm = 1.0;
p.y_O2 = 0.21;
p.y_CO2_atm = 3.8e-4;
p.y_CO2_inj = 1.0;
p.KH_O2_ref = 1.3e-3 * 1e3;
p.KH_CO2_ref = 3.4e-2 * 1e3;
p.C_O2 = 1700;
p.C_CO2 = 2400;
p.K1_ref = 4.3e-7 * 1e3;
p.K2_ref = 4.7e-11 * 1e3;
p.KW_ref = 1e-14 * 1e6;
p.dH_K1 = -7500;
p.dH_K2 = -15000;
p.dH_KW = 55900;

p.alpha_CO2 = 0.05;
p.beta_CO2 = 0.7;
p.alpha_O2 = 0.04;
p.beta_O2 = 0.7;
p.k_strip_CO2_O2 = 0.05;
p.k_strip_O2_CO2 = 0.02;
p.k_atm_CO2 = 1e-5;
p.k_atm_O2 = 2e-5;
p.k_pw_CO2 = 1e-5;
p.k_pw_O2 = 2e-5;
p.DIC_in = 2;
p.Cat_in = 5;

end

function [mu_I, mu_T, mu_pH, mu_DO, mu, m, P] = i_growth_terms(Xalg, XO2, ~, H, T, Depth, Env, p)
N = numel(Xalg);
PAR = 0.46 * 4.56 * Env.RadG;

Iav = zeros(N,1);
for k = 1:N
    den = p.Ka * max(Depth(k), 1e-6) * max(Xalg(k), 1e-8);
    Iav(k) = PAR(k) / max(den, 1e-8) * (1 - exp(-den));
end

mu_I = Iav.^p.n_hill ./ (p.Ik^p.n_hill + Iav.^p.n_hill + 1e-12);
mu_T = i_window_numeric(T, p.T_min, p.T_opt, p.T_max);
pH = -log10(max(H, 1e-12) / 1000);
mu_pH = i_window_numeric(pH, p.pH_min, p.pH_opt, p.pH_max);

KH_O2 = p.KH_O2_ref * exp(p.C_O2 * (1./(T+273.15) - 1/p.T_ref));
Xeq_O2 = KH_O2 * p.p_atm * p.y_O2;
DO_sat = 100 * XO2 ./ max(Xeq_O2, 1e-12);
mu_DO = max(0, 1 - (DO_sat / p.DO_max).^p.m_DO);

P = p.mu_max .* mu_I .* mu_T .* mu_pH .* mu_DO;
mu = p.eta_X .* P;
m = p.m_min .* (1 + p.k_resp_I .* (1 - mu_I)) .* p.Q10.^((T - 20) / 10);

end

function w = i_window_numeric(x, a, b, c)
w = zeros(size(x));

idx1 = x > a & x <= b;
r1 = (x(idx1) - a) / max(b - a, eps);
w(idx1) = 3*r1.^2 - 2*r1.^3;

idx2 = x > b & x < c;
r2 = (c - x(idx2)) / max(c - b, eps);
w(idx2) = 3*r2.^2 - 2*r2.^3;

w = min(max(w, 0), 1);

end

function [CO2, HCO3, CO3] = i_carbonate_species(DIC, H, T, p)
K1 = p.K1_ref * exp(-p.dH_K1/p.R_gas * (1./(T+273.15) - 1/p.T_ref));
K2 = p.K2_ref * exp(-p.dH_K2/p.R_gas * (1./(T+273.15) - 1/p.T_ref));

Delta = H.^2 + H.*K1 + K1.*K2;
CO2 = DIC .* H.^2 ./ max(Delta, 1e-16);
HCO3 = DIC .* H .* K1 ./ max(Delta, 1e-16);
CO3 = DIC .* K1 .* K2 ./ max(Delta, 1e-16);

end

function [t, Env] = i_build_time_and_env(Data)
RadG = [];
RadPAR = [];
Temp_ext = [];
RH = [];
Wind = [];

for i = 1:numel(Data)
    if ~isfield(Data(i), 'u') || isempty(Data(i).u)
        continue;
    end

    ui = double(Data(i).u);
    if size(ui,2) < 5 && size(ui,1) >= 5
        ui = ui.';
    end
    if size(ui,2) < 5
        error('Data(%d).u must have at least 5 columns.', i);
    end

    RadG = [RadG; ui(:,1)]; %#ok<AGROW>
    RadPAR = [RadPAR; ui(:,2)]; %#ok<AGROW>
    Temp_ext = [Temp_ext; ui(:,3)]; %#ok<AGROW>
    RH = [RH; ui(:,4)]; %#ok<AGROW>
    Wind = [Wind; ui(:,5)]; %#ok<AGROW>
end

N = numel(RadG);
t = (0:N-1).' * 60;

Env = struct();
Env.RadG = RadG;
Env.RadPAR = RadPAR;
Env.Temp_ext = Temp_ext;
Env.RH = RH;
Env.Wind = Wind;

end

function future = i_build_future(t, Env, k)
h = min(60, numel(t) - k);
idx = (k+1):(k+h);
if isempty(idx)
    idx = k;
end

future = struct();
future.t_future = t(idx);
future.RadGlobal = Env.RadG(idx);
future.RadPAR = Env.RadPAR(idx);
future.Temp_ext = Env.Temp_ext(idx);
future.RH = Env.RH(idx);
future.Wind = Env.Wind(idx);

end

function v = i_getfield_or(s, fname, default_value)
if isstruct(s) && isfield(s, fname) && ~isempty(s.(fname))
    v = s.(fname);
else
    v = default_value;
end
end
