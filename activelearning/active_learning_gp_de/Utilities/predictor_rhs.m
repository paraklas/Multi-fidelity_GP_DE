function [f, v] = predictor_rhs(x_star)

global ModelInfo

x0 = ModelInfo.x0;
x1 = ModelInfo.x1;

u0 = ModelInfo.u0;
y1 = ModelInfo.y1;

hyp = ModelInfo.hyp;

D = size(x1,2);

y = [u0; y1];

L=ModelInfo.L;

psi1 = h(x0, x_star, hyp(1:D+1),0)';
psi2 = k(x1, x_star, hyp(1:D+1),0)';
psi = [psi1 psi2];

% calculate prediction
f = psi*(L'\(L\y));

v = k(x_star, x_star, hyp(1:D+1),0) ...
  - psi*(L'\(L\psi'));

v = abs(diag(v));