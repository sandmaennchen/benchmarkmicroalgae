function sim_solver = get_acados_simsolver

model = get_acados_raceway_model();

sim = AcadosSim();
sim.model = model;
sim.solver_options.Tsim = 60;
sim.solver_options.integrator_type = 'ERK';


sim_solver = AcadosSimSolver(sim);


end