function [f, v] = predictor(x_star)

global ModelInfo

x0 = ModelInfo.x0;
x1 = ModelInfo.x1;

u0 = ModelInfo.u0;
y1 = ModelInfo.y1;

hyp = ModelInfo.hyp;

D = size(x1,2);

y = [u0; y1];

L=ModelInfo.L;

psi1 = g(x_star, x0, hyp(1:D+1),0);
psi2 = h(x_star, x1, hyp(1:D+1),0);
psi = [psi1 psi2];

% calculate prediction
f = psi*(L'\(L\y));

v = g(x_star, x_star, hyp(1:D+1),0) ...
  - psi*(L'\(L\psi'));

v = abs(diag(v));