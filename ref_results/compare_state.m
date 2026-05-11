data_ref = load('state_traj_ref.mat');
traj_ref = data_ref.state_traj;
data = load('state_traj.mat');
traj = data.state_traj;

N = length(traj_ref);
for n = 1:N
    val_ref = traj_ref(n);
    val = traj(n);

    if isequal(val, val_ref)
        fprintf('state is identical at step %d\n', n);
    else
        fprintf('state differs at step %d\n', n);

        val_ref{1}
        val{1}

        keyboard

    end
end