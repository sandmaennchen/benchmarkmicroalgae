function [results] = simulate_benchmark_model_approximation(Data, ctrl)
% Approximate benchmark simulator with compatible interface.
% The plant update is intentionally simplified (identity mapping) so users
% can test controller wiring and result visualization without the black-box model.

if nargin < 2
    error('simulate_benchmark_model_approximation requires Data and ctrl inputs.');
end

% Build simulation horizon and environmental series from Data.
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

refs = struct('pH', 8.0, 'DO', 150.0, 'T', 25.0);

% Identity-mapped states (constant over time in this approximation).
pH0    = refs.pH;
DO0    = refs.DO;
Depth0 = 0.15;
X0_gL  = 0.50;
T0     = Env.Temp_ext(1);

pH    = pH0    * ones(N,1);
DO    = DO0    * ones(N,1);
Depth = Depth0 * ones(N,1);
X_gL  = X0_gL  * ones(N,1);
T     = T0     * ones(N,1);

QCO2_cmd = zeros(N,1);
Qair_cmd = zeros(N,1);
QCO2_del = zeros(N,1);
Qair_del = zeros(N,1);
Qd_bin   = zeros(N,1);
Qh_bin   = zeros(N,1);
Qd       = zeros(N,1);
Qh       = zeros(N,1);
Qhx_cmd  = zeros(N,1);
Tin_hx   = T0 * ones(N,1);

st_pH = struct();
st_DO = struct();
st_HD = struct();
st_TX = struct();

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
        [st_cmd, st_pH] = ctrl.fn_pH_CO2(Timeline, obs, refs, env, future, st_cmd, st_pH);
    end
    if isfield(ctrl, 'fn_DO_air') && ~isempty(ctrl.fn_DO_air)
        [st_cmd, st_DO] = ctrl.fn_DO_air(Timeline, obs, refs, env, future, st_cmd, st_DO);
    end
    if isfield(ctrl, 'fn_HD') && ~isempty(ctrl.fn_HD)
        [st_cmd, st_HD] = ctrl.fn_HD(Timeline, obs, refs, env, future, st_cmd, st_HD);
    end
    if isfield(ctrl, 'fn_Temp_HX') && ~isempty(ctrl.fn_Temp_HX)
        [st_cmd, st_TX] = ctrl.fn_Temp_HX(Timeline, obs, refs, env, future, st_cmd, st_TX);
    end

    QCO2_cmd(k) = i_getfield_or(st_cmd, 'Qco2', 0);
    Qair_cmd(k) = i_getfield_or(st_cmd, 'Qair', 0);
    Qd_bin(k)   = i_getfield_or(st_cmd, 'Qd_bin', 0);
    Qh_bin(k)   = i_getfield_or(st_cmd, 'Qh_bin', 0);
    Qhx_cmd(k)  = i_getfield_or(st_cmd, 'Qhx', 0);
    Tin_hx(k)   = i_getfield_or(st_cmd, 'Tin_hx', obs.T);
end

% Simple actuator limits and binary-flow conversion.
QCO2_max = 20/1000/60;   % [m^3/s], 20 L/min
Qair_max = 500/1000/60;  % [m^3/s], 500 L/min
Qd_rate  = 20/1000/60;   % [m^3/s], approximate fixed dilution flow
Qh_rate  = 20/1000/60;   % [m^3/s], approximate fixed harvest flow

QCO2_del = min(max(QCO2_cmd, 0), QCO2_max);
Qair_del = min(max(Qair_cmd, 0), Qair_max);
Qd       = Qd_rate * (Qd_bin > 0.5);
Qh       = Qh_rate * (Qh_bin > 0.5);

% Cumulative usage.
cum_air_L = cumsum(Qair_del * dt * 1000);
cum_CO2_L = cumsum(QCO2_del * dt * 1000);

% Coarse KPI placeholders (consistent dimensions/units for plotting/reporting).
area_m2 = 1.0;
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

% Limitation/growth placeholders.
mu_I  = Env.RadPAR / max(max(Env.RadPAR), 1);
mu_T  = ones(N,1);
mu_pH = ones(N,1);
mu_DO = ones(N,1);
mu    = zeros(N,1);
m     = zeros(N,1);
P     = zeros(N,1);

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

results.DIC = zeros(N,1);
results.Cat = zeros(N,1);
results.HCO3 = zeros(N,1);
results.CO3 = zeros(N,1);
results.CO2 = zeros(N,1);

results.HX = struct();
results.HX.Qw_m3s = max(Qhx_cmd, 0);
results.HX.Tin_C = Tin_hx;
results.HX.Tout_C = Tin_hx;
results.HX.UA_WK = zeros(N,1);
results.HX.Q_W = zeros(N,1);
results.HX.limits = struct('Qw_min',0,'Qw_max',1e-2,'Tin_min',0,'Tin_max',50);

results.Env = Env;

results.J = struct();
results.J.pH = mean((pH - refs.pH).^2);
results.J.DO = mean((DO - refs.DO).^2);
results.J.T  = mean((T - refs.T).^2);

end

function [t, Env] = i_build_time_and_env(Data)
% Assemble a continuous minute-wise horizon from the daily data blocks.

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
% Provide a short look-ahead trajectory for controllers that can use it.

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