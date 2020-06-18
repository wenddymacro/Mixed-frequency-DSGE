var y, r, pii // endogenous variables
z, v;        // shock processes

varexo eps_z, eps_v; // exogenous shocks

parameters rho_r, phi_pii, phi_y, phii,  thetta, betta // model parameters
rho_z, rho_v                                        // shock parameters
ypot, kappa;            

// parameter values
betta = 0.96;  
phi_y = 1.5;    
phi_pii = 8;
phii = 1;
rho_r = 0.8;
thetta = 0.8;
ypot = 0;

// shock processes
rho_z = 0.9;
rho_v = 0.4;

kappa = ((1-thetta)*(1-betta*thetta)/thetta)*(1+phii);


model(linear);

// structural equations
y = y(+1) - (r - pii(+1)+ log(betta)  ) + (1-rho_z)*z; // IS-Curve
r = rho_r*r(-1) + phi_pii*pii(+1) + phi_y*(y-ypot) + v;    // Taylor Rule
pii = betta*pii(+1) + ((1-thetta)*(1-betta*thetta)/thetta)*(1+phii) * (y-ypot); // NK-PC

// shock processes
z = rho_z*z(-1) + eps_z;
v = rho_v*v(-1) + eps_v;
end;

varobs y r;
check;
resid;


shocks;
var eps_z; stderr 0.01;
var eps_v; stderr 0.01;
end;

estimated_params;
betta, 0.96 ,0.9, 1.0;          % beta
thetta, 0.8 ,0.1, 0.9;         % kappa
rho_r, 0.8 ,0.5, 0.99;            % rho
phi_y, 1.5 ,0, 3;              % phiy
phi_pii, 8 ,5, 10;
rho_v, 0.4 ,0.1, 0.99;
rho_z, 0.9 ,0.5, 0.99;
end;

estimation(datafile = Toy_data, xls_sheet = Sheet1, mh_replic = 40000);

//stoch_simul(irf=40);
//stoch_simul(periods=2000);
