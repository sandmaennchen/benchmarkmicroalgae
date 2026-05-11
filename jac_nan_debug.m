clear;
close all;
clc;

import casadi.*

% Symbolics
x = SX.sym('x', 7);   % [Xalg, XO2, DIC, Cat, H, T, V]
u = SX.sym('u', 6);   % [QCO2, Qair, Qd, Qh, Qw, Tin_hx]
d = SX.sym('d', 4);   % [RadG, RH, Text, Uwind]

% Some constants
QCO2_max = 20/1000/60;
Qair_max = 500/1000/60;
Qd_rate  = 20/1000/60;
Qh_rate  = 20/1000/60;
Qw_min = 0;
Qw_max = 1e-2;
Tin_min = 0;
Tin_max = 50;

% Data and parammeters
load_data;
p = create_p();

% Build integrator
dt = 60; % time step
[ode_fun, out_fun] = build_ode(p);
dae = struct('x', x, 'p', vertcat(u, d), 'ode', ode_fun(x, u, d));
opts = struct('tf', dt, 'abstol', 1e-5, 'reltol', 1e-6);
F_step = integrator('F_step', 'cvodes', dae, opts);

% Initial states
[t, Env] = build_time_and_env(Data);
x0 = get_initial_state(p, Env);
u0 = [QCO2_max; Qair_max; Qd_rate; Qh_rate; Qw_max; Tin_max]; % !Guessed
d0 = [Env.RadG(1); Env.RH(1); Env.Temp_ext(1); Env.Wind(1)];

% Define jacobians functions
ode_fun_jac = Function('ode_fun_jac', {x, u, d} , {jacobian(ode_fun(x, u, d), x)}); % Jacobian of ODE_FUNCTION
% ode_symb_jac = Function('ode_jac', {x, u, d} , {jacobian(xdot, x)}); % Jacobian of ODE_SYMBOLIC

% Debug
% ode_fun_jac(x, u, d)
% ode_fun_jac(x0, u0, d0)
% 
% ode_symb_jac(x, u, d)
% ode_symb_jac(x0, u0, d0)

out = F_step('x0', x0, 'p', [u0; d0]);
out.xf;

N = 1440;
p_input = [u0; d0];
p_batch = [QCO2_max; Qair_max; 0; 0; 0; Tin_max; d0];
% p_problamatic = [0; 0; 0; 0; 0; 38.7125; 0.28125; 24.5625; 38.7125; 0.0533796]; % This causes error when using collocations, works with cvodes
P_seq = repmat(p_batch, 1, N);
F_step10 = F_step.mapaccum('F_step10', N);
out10 = F_step10('x0', x0, 'p', P_seq);
out10.xf;
% plot(full(out10.xf(6,:)))
plot(-log10(full(out10.xf(5,:))/1000))