function x0 = get_initial_state(p, Env)
Xalg0 = 0.5 * 1e3; % tv: ref results 0.5; TABLE: 0.32 [g/l]


% From eq 4 eq 57 and eq 60 (DO = 89)
XO20 = p.KH_O2_ref * p.p_atm * p.y_O2;


% DIC0 = 0.014 * 1e3; % v
% Cat0 = 1.672669964353326e+01; % p.Cat_in;
DIC0 = 17.5;
Cat0 = 17.0;
H0 = 1.582e-05; % v
T0 = 25;
Depth0 = 0.15;
V0 = max(p.Vsump + 1e-2, p.A * Depth0 + p.Vsump);
x0 = [Xalg0; XO20; DIC0; Cat0; H0; T0; V0];

end