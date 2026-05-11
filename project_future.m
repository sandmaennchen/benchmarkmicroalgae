function future = project_future(t, Env, k)
h = min(60, numel(t) - k);
idx = (k+1):(k+h);
if isempty(idx)
    idx = k;
end

future = struct();
future.t_future = t(idx);
future.RadGlobal = Env.RadG(idx);
future.RadPAR = Env.RadPAR(idx);
future.Temp_ext = Env.Temp_ext(idx);
future.RH = Env.RH(idx);
future.Wind = Env.Wind(idx);

end