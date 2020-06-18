var y, pii, R, z, nu;

varexo eps_z, eps_nu;

parameters kapa, phi_y, rho_r, betta, phi_pi, rho_nu, rho_z, y_pot;

betta = 0.96;          % beta
kapa = 0.0907;         % kappa
rho_r = 0.85;            % rho
phi_y = 1;              % phiy
phi_pi = 7.6359;
rho_nu = 0.8;
rho_z = 0.8;
y_pot = 0;

model(linear);

y = y(+1) - (R - pii(+1) + log(betta)) + (1 - rho_z) * z;
R = rho_r * R(-1) + phi_pi * pii(+1) + phi_y * (y - y_pot) + nu;
pii = betta * pii(+1) + kapa * (y - y_pot);

nu = rho_nu * nu(-1) + eps_nu;
z = rho_z * z(-1) + eps_z;

end;


varobs y R;
check;
resid;

initval;
y = 0;
pii = 1;
R = 1;
end;

shocks;
var eps_nu; stderr 0.01;
var eps_z; stderr 0.01;
end;



estimated_params;
betta, 0.96 ,0.9, 1.0;          % beta
kapa, 0.0907 ,0, 0.5;         % kappa
rho_r, 0.85 ,0.5, 0.99;            % rho
phi_y, 1 ,0, 3;              % phiy
phi_pi, 7.6359 ,5.1, 10;
rho_nu, 0.8 ,0.5, 0.99;
rho_z, 0.8 ,0.5, 0.99;
end;


estimation(datafile = Toy_data, xls_sheet = Sheet1, mh_replic = 40000);

