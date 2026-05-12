function xn = acados_integrator(sim_solver, xk, uk, dk)

sim_solver.set('x', xk);
sim_solver.set('u', uk);
sim_solver.set('p', dk);
sim_solver.solve();

xn = sim_solver.get('xn');

end