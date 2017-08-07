function [f, v] = predictor(x_star)

global ModelInfo

x0 = ModelInfo.x0;
x1 = ModelInfo.x1;
x2 = ModelInfo.x2;
u0 = ModelInfo.u0;
y1 = ModelInfo.y1;
y2 = ModelInfo.y2;
hyp = ModelInfo.hyp;
rho = hyp(end-3);

D = size(x2,2);

y = [u0; y1; y2];

L=ModelInfo.L;

psi1 = rho^2*g(x_star, x0, hyp(1:D+1),0) ...
           + g(x_star, x0, hyp(D+2:2*D+2),0);
psi2 = rho*h(x_star, x1, hyp(1:D+1),0);
psi3 = rho^2*h(x_star, x2, hyp(1:D+1),0) + h(x_star, x2, hyp(D+2:2*D+2),0);
psi = [psi1 psi2 psi3];

% calculate prediction
f = psi*(L'\(L\y));

v = rho^2*g(x_star, x_star, hyp(1:D+1),0) ...
  + g(x_star, x_star, hyp(D+2:2*D+2),0) ...
  - psi*(L'\(L\psi'));

v = diag(v);