function F_step = build_integrator(ode_fun, dt)
import casadi.*

x = SX.sym('x', 7);
u = SX.sym('u', 6);
d = SX.sym('d', 4);
xdot = ode_fun(x, u, d);

dae = struct('x', x, 'p', vertcat(u, d), 'ode', xdot);
opts = struct('tf', dt, 'abstol', 1e-5, 'reltol', 1e-6);
F_step = integrator('F_step', 'cvodes', dae, opts);
% F_step = integrator('F_step', 'collocation', dae, 0, dt, struct('collocation_scheme', 'radau', 'interpolation_order', 2, 'rootfinder', 'kinsol'));

end