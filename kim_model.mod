var a, i, k, y, l, c, w, r;

varexo eps_a, eps_y, eps_c, eps_i, eps_r, eps_l, eps_w;

parameters delta, alpha, CY, YK, sigma_n, sigma_c, beta, rho_a, R, sigma_z, K, L;

delta = 0.025/3;
alpha = 0.3;
sigma_n = -1; 
sigma_c = 1; 
beta = 0.98;
rho_a = 0.9;
R = 1/beta - 1 + delta;
sigma_z = 1;
L = 0.2;
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

end;

shocks;
var eps_a; stderr sigma_z;
var eps_y; stderr 1.3964;
var eps_c; stderr 0.9529;
var eps_w; stderr 1.0031;
var eps_l; stderr 0.6637;
var eps_r; stderr 1.1439;
var eps_i; stderr 9.7956;
end;


steady;
check; 
resid;

stoch_simul(hp_filter = 14400, nograph, periods = 300) r w l y c i;