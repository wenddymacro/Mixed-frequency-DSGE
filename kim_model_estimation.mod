var i, y, c, w, r, l, a, k;

varexo eps_a, eps_y, eps_c, eps_i, eps_r, eps_l, eps_w;

parameters delta, alpha, CY, YK, sigma_n, sigma_c, beta, rho_a, R, sigma_z, K, L;

delta = 0.025/3;
alpha = 0.3;
sigma_n = -1.00001; 
sigma_c = 1; 
beta = 0.98;
rho_a = 0.9;
sigma_z = 1;
L = 0.2;
R = 1/beta - 1 + delta;
K = ( R / (alpha * L ^ (1-alpha) ) ) ^ ( 1 / (alpha - 1) );
CY = 1 - delta * (K/L) ^ (1-alpha);
YK = delta / (1-CY);

model(linear);

a = rho_a * a(-1) + eps_a;
i = 1/delta * (k - (1-delta) * k(-1)) + eps_i;
y = a + alpha * k(-1) + (1-alpha) * l + eps_y;
k = YK * y - YK * CY * c + (1-delta) * k(-1);
w = a + alpha * k(-1) - alpha * l + eps_w;
r = a + (alpha - 1) * k(-1) + (1-alpha) * l + eps_r;
l = 1/sigma_n * (-sigma_c * c + w) + eps_l;
c = -beta * R / sigma_c * r(+1) + c(+1) + eps_c;
//c = beta * R / sigma_c * r + c(-1) + eps_c;
end;


varobs i, y, c, w, r, l;
//check;
//resid;

steady;
check;
stoch_simul;

shocks;
var eps_a; stderr sigma_z;
var eps_y; stderr 1.3964;
var eps_c; stderr 0.9529;
var eps_w; stderr 1.0031;
var eps_l; stderr 0.6637;
var eps_r; stderr 1.1439;
var eps_i; stderr 9.7956;
end;


estimated_params;
alpha, uniform_pdf,,,0.2, 0.4;         
beta, uniform_pdf,,,0.97, 0.99;            
delta, uniform_pdf,,,0.006, 0.010;              
rho_a, uniform_pdf,,,0.8, 1;
stderr eps_a, inv_gamma_pdf, 1, 0.12;
stderr eps_r, inv_gamma_pdf, 1.2, 0.12;
stderr eps_l, inv_gamma_pdf, 0.6, 0.05;
stderr eps_w, inv_gamma_pdf, 1, 0.12;
stderr eps_y, inv_gamma_pdf, 1.5, 0.25;
stderr eps_c, inv_gamma_pdf, 1, 0.12;
stderr eps_i, inv_gamma_pdf, 10, 2;
end;

estimated_params_init(use_calibration);
end;



estimation(datafile = Kim_data_mixed, xls_sheet = Sheet1, mh_replic = 2000000, mh_drop = 0.5) y r c i a l;
//stoch_simul(hp_filter = 14400, nograph, periods = 300);