function ocp_solver = get_acados_algae_ocp_solver()

p = create_p()
ocp = AcadosOcp();
model = get_acados_raceway_model()
ocp.model = model;

% Notes OCP:
% - model;
%  - works. but might need modifications to be well scaled
%  - H bad scaled -> use pH instead
%  - Cat -> pOH
%  - rates of controls needed. not clear if for all.

% - objective:
%  J_SP + J_SS + J_RC
%  setpoint tracking; control signal smoothness (rate cost); resource consumptoin

% - constraints: physical ones, control bounds.



% States
Xalg = model.x(1);
XO2  = model.x(2);
DIC  = model.x(3);
Cat  = model.x(4);
H    = model.x(5);
T    = model.x(6);
V    = model.x(7);

% Inputs
QCO2   = model.u(1);
Qair   = model.u(2);
Qd     = model.u(3);
Qh     = model.u(4); % harvest.
Qw     = model.u(5);
Tin_hx = model.u(6);  % Heat ex. command
% -> this goes to plant


T_K = T + 273.15;
KH_O2  = p.KH_O2_ref  * exp(p.C_O2  * (1/T_K - 1/p.T_ref)); % equation (60)
% KH_CO2 = p.KH_CO2_ref * exp(p.C_CO2 * (1/T_K - 1/p.T_ref)); % equation (60)
Xeq_O2   = KH_O2 * p.p_atm * p.y_O2; % equation (57)

% TODO:
pH = -log(H / 1000) / log(10); % equation (3)
DO = 100 * XO2 / (Xeq_O2 + p.reg_eps); % equation (4)
Xalg_gL = Xalg / 1000; % equation (5)
depth = (V - p.Vsump) / (p.W * p.L + p.reg_eps); % equation (6)

% Control bounds
Qd_rate  = 20/1000/60;
Qh_rate  = 20/1000/60;
QCO2_rate = 20/1000/60;
Qair_rate = 500/1000/60;

QCO2_min = 0
Qair_min = 0
Qd_min = 0;
Qh_min = 0;

Qd_max = 1;
Qh_max = 1;
QCO2_max = 1;
Qair_max = 1;

Qw_min = 0;
Qw_max = 1e-2;
Tin_min = 0; % TODO: might be integer?!
Tin_max = 50;

nu = length(model.u);



ocp.constraints.idxbu = 0:nu-1
ocp.constraints.lbu = [QCO2_min,
                       Qair_min,
                       Qd_min,
                       Qh_min,
                       Qw_min,
                       Tin_min]
ocp.constraints.ubu = [QCO2_max,
                       Qair_max,
                       Qd_max,
                       Qh_max,
                       Qw_max,
                       Tin_max]

% intial state constraint // BIT DIRTY
load_data;
[t, Env] = build_time_and_env(Data);
x0 = get_initial_state(p, Env);
ocp.constraints.x0 = x0;


% Cost
pH_ref = 8.0;
dO_ref = 150; % [%]
T_ref = 30; % [C]
depth_ref = .15; % [m]

% define nonlinear cost, later convex-over-nonlinear might be better for convergence.
cost = 0;

% setpoint tracking; NOTE: paper defines l1-tracking normalized by ref value
cost = cost + 1 * (pH - pH_ref)^2;
cost = cost + 1 * (DO - dO_ref)^2;
cost = cost + 1 * (T - T_ref)^2;
% cost = cost + 0.01 * (depth - depth_ref)^2;

% % control effort
% cost = cost + 1 * QCO2 + 1 * Qair + 1 * Qd + 1 * Qh + 1 * Qw + 1 * Tin_hx;
% % TODO: smoothness cost; augment rate.

% % harvested amount:
% cost = cost + -1 * (Xalg_gL * depth * p.A) * Qh_rate * Qh; % algae amount * harvest rate

% cost_e = -1 * Xalg_gL * depth * p.A; % remaining amount of algae.
cost_e = 0;

ocp.model.cost_expr_ext_cost = cost;
ocp.model.cost_expr_ext_cost_e = cost_e;

ocp.solver_options.tf = 60*60*24;
dt = 60*5;
N_horizon = ocp.solver_options.tf/dt;
ocp.solver_options.N_horizon = N_horizon;
ocp.solver_options.integrator_type = 'ERK';
ocp.solver_options.hessian_approx = 'EXACT';
ocp.solver_options.regularize_method = 'MIRROR';
ocp.solver_options.reg_epsilon = 1e-2;
ocp.solver_options.globalization = 'MERIT_BACKTRACKING';

ocp.cost.cost_type = 'EXTERNAL';
ocp.cost.cost_type_e = 'EXTERNAL';


ocp_solver = AcadosOcpSolver(ocp);

% test call
ocp_solver.set('constr_x0', x0);


% parameters:
i_cl = 0;
for k=0:N_horizon
    i_env = i_cl + 5*k + 1;
    dk = [Env.RadG(i_env); Env.RH(i_env); Env.Temp_ext(i_env); Env.Wind(i_env)];
    ocp_solver.set('p', dk, k);
end

% initialize u;
for k=0:N_horizon-1
    u_guess = 1e-6 * ones(nu, 1);
    u_guess(nu) = 25; %Env.Temp_ext(1);
    ocp_solver.set('u', u_guess, k);
end

ocp_solver.solve();
ocp_solver.print_statistics();
status = ocp_solver.get('status')
ocp_solver.get('qp_Q', 1)
ocp_solver.get('qp_R', 1)
ocp_solver.dump_last_qp_to_json('Matlab')

keyboard

% if successfully solved, analyze solutions
ocp_solver.get('u')
ocp_solver.get('x')


end