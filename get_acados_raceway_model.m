function model = get_acados_raceway_model()
import casadi.*

x = SX.sym('x', 7);
u = SX.sym('u', 6);
d = SX.sym('d', 4);
xdot = SX.sym('xdot', 7);

p = create_p();
[ode_fun, ~, ~] = build_ode(p);

f_expl = ode_fun(x, u, d);

model = AcadosModel();
model.name = 'raceway_model';
model.xdot = xdot;
model.x = x;
model.u = u;
model.p = d;
model.f_expl_expr = f_expl;
model.f_impl_expr = f_expl - xdot;

end