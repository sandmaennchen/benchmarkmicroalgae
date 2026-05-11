function [t, Env] = build_time_and_env(Data)
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
    Temp_ext = [Temp_ext; ui(:,4)]; %#ok<AGROW>
    RH = [RH; ui(:,3)]; %#ok<AGROW>
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