var y, pii, R;

varexo eps_rt, eps_dt, eps_st;

parameters kapa, phi_y, rho_r, betta, phi_pi, tau, theta;

betta = 0.99;           % beta
rho_r = 0.9;            % rho
phi_y = 0.5;              % phiy
phi_pi = 1.5;
tau = 1;
theta = 0.9;

kapa = (1-betta*theta) * (1-theta) / theta;

model(linear);

y = y(+1) - tau * (R - pii(+1)) + eps_dt;
R = rho_r * R(-1) + (1-rho_r) * (phi_pi * pii + phi_y * y) + eps_rt;
pii = betta * pii(+1) + kapa * y + eps_st;



end;



varobs y pii R;
check;
resid;

estimated_params;
theta, 0.9 ,0.5, 1;         % kappa
rho_r, 0.9 ,0.5, 0.99;            % rho
phi_y, 0.5 ,0, 3;              % phiy
phi_pi, 1.5 ,1, 3;
stderr eps_st, 0.1, 0.01, 1;
stderr eps_dt, 0.1, 0.01, 1;
stderr eps_rt,  0.1, 0.01, 1;
end;



estimation(datafile = Foroni_data_mixed, xls_sheet = Sheet1, mh_replic = 1000) y pii R;
stoch_simul(hp_filter = 14400, nograph, periods = 300);