var y, pii, R;

varexo eps_rt, eps_dt, eps_st;

parameters kapa, phi_y, rho_r, betta, phi_pi, tau, theta;

betta = 0.99;           % beta
//rho_r = 0.9;            % rho
//phi_y = 0.5;              % phiy
//phi_pi = 1.5;
tau = 1;
//theta = 0.9;

theta = 0.9;
rho_r = 0.5;
phi_y = 1.3069;
phi_pi = 1.7739;

kapa = (1-betta*theta) * (1-theta) / theta;

model(linear);

y = y(+1) - tau * (R - pii(+1)) + eps_dt;
R = rho_r * R(-1) + (1-rho_r) * (phi_pi * pii + phi_y * y) + eps_rt;
pii = betta * pii(+1) + kapa * y + eps_st;


end;


shocks;
//var eps_dt; stderr 0.1;
//var eps_rt; stderr 0.1;
//var eps_st; stderr 0.1;

var eps_dt; stderr 1;
var eps_rt; stderr 0.6385;
var eps_st; stderr 0.3539;
end;

steady;
check; 
resid;
stoch_simul(hp_filter = 14400, nograph, periods = 300) y pii R;
