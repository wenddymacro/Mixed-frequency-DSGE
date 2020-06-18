var y, g, r, pii, z, dyobs, robs, pinfobs;

varexo eps_r, eps_z, eps_g;

parameters tau, rho_z, betta, kappa, psi1, psi2, rho_r, rho_g, gamma, pistar, rstar, sigma_z, sigma_g, sigma_r;

betta = 0.99; 
gamma = 0.5;
pistar = 5;
rstar = 2;
tau = 2;
kappa = 0.3;
psi1 = 1.5;
psi2 = 0.125;
rho_g = 0.9;
rho_z = 0.2;
rho_r = 0.7;
sigma_g = 1.25;
sigma_z = 1.25;
sigma_r = 0.63;

model(linear);

y = g + (y(+1) - g(+1)) - 1 / tau * (r - pii(+1) - rho_z * z);
pii = betta * pii(+1) + kappa * (y - g);
r = psi1 * (1 - rho_r) * pii + psi2 * (1 - rho_r) * y + rho_r * r(-1) + sigma_r * eps_r;
z = rho_z * z(-1) + sigma_z * eps_z;
g = rho_g * g(-1) + sigma_g * eps_g;
dyobs = ln(gamma) + y - y(-1) + z;
robs = pistar + rstar + 4 * r;
pinfobs = pistar + 4 * pii;
end;


varobs dyobs, robs, pinfobs;
check;
resid(1);

steady;
check;
stoch_simul;

shocks;
var eps_r; stderr 1;
var eps_z; stderr 1;
var eps_g; stderr 1;
end;


estimated_params;
gamma, normal_pdf, 0.5, 0.5;
pistar, gamma_pdf, 5, 2;
rstar, gamma_pdf, 2, 1;
tau, gamma_pdf, 2, 0.5;
kappa, gamma_pdf, 0.3, 0.1;
psi1, gamma_pdf, 1.5, 0.5;
psi2, gamma_pdf, 0.125, 0.1;
rho_g, beta_pdf, 0.9, 0.05;
rho_z, beta_pdf, 0.2, 0.1;
rho_r, beta_pdf, 0.7, 0.15;
sigma_g, inv_gamma_pdf, 1.25, 0.65;
sigma_z, inv_gamma_pdf, 1.25, 0.65;
sigma_r, inv_gamma_pdf, 0.63, 0.33;
end;

estimated_params_init(use_calibration);
end;

estimation(datafile = Monti_replication, xls_sheet = Sheet1, mh_replic = 2000000, mh_drop = 0.5);
//stoch_simul(hp_filter = 14400, nograph, periods = 300);