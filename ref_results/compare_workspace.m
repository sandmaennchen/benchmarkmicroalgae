workspace_ref = load('ref_results.mat');
workspace = load('matlab.mat');

vars = fieldnames(workspace_ref);

for i = 1:numel(vars)
    v = vars{i};
    
    if isequal(workspace_ref.(v), workspace.(v))
        fprintf('Variable %s is identical\n', v);
    else
        fprintf('Variable %s differs\n', v);
    end
end

control_vars = fieldnames(workspace_ref.ctrl);

for i = 1:numel(control_vars)
    v = control_vars{i};
    
    if isequal(workspace_ref.ctrl.(v), workspace.ctrl.(v))
        fprintf('Variable %s is identical\n', v);
    else
        fprintf('Variable %s differs\n', v);
    end
end