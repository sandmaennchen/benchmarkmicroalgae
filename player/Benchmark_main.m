%% ===== Benchmark Microalgae Raceway Reactor =====
clear; 
close all;
clc;
rng(42,'twister');
addpath('..');

%% ---------------- Data load ----------------
load_data;


global obs_traj
obs_traj = {};
global state_traj
state_traj = {};
%% ---------------- Controllers (one per actuator) ----------------
ctrl = struct();
ctrl.fn_pH_CO2  = @controller_pH_OnOff;
ctrl.fn_DO_air  = @controller_DO_OnOff;
ctrl.fn_HD      = @controller_HD_fixed;
ctrl.fn_Temp_HX = @controller_Temp_HX_no_control;

%% ---------------- Simulation ----------------
fprintf('\nRunning Benchmark Microalgae Raceway simulation...\n');
% [results] = simulate_benchmark_model(Data, ctrl);
[results] = simulate_benchmark_model_approximation(Data, ctrl);

save('state_traj_ref.mat', 'state_traj')
save('obs_traj_ref.mat', 'obs_traj')


%% ---------------- Results ----------------
show_results;