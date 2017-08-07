function [f, v] = predictor_rhs_H(x_star)

global ModelInfo
hyp = ModelInfo.hyp;
x1 = ModelInfo.x1;
x2 = ModelInfo.x2;

y1 = ModelInfo.y1;
y2 = ModelInfo.y2;
y=[y1;y2];

jitter = ModelInfo.jitter;

sigma_n2 = exp(hyp(end));
sigma_n1 = exp(hyp(end-1));
rho = hyp(end-3);

n1 = size(x1,1);
[n2,D] = size(x2);
n = n1+n2;

K11 = k(x1, x1, hyp(1:D+1),0);
K12 = rho*k(x1, x2, hyp(1:D+1),0);
K21 = rho*k(x2, x1, hyp(1:D+1),0);
K22 = (rho^2)*k(x2, x2, hyp(1:D+1),0) + k(x2, x2, hyp(D+2:2*D+2),0);

K11 = K11 + eye(n1).*sigma_n1;
K22 = K22 + eye(n2).*sigma_n2;

K = [K11 K12;
    K21 K22];

K = K + eye(n).*jitter;

% Cholesky factorisation
[L,p]=chol(K,'lower');

if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

psi1 = rho*k(x_star,x1,ModelInfo.hyp(1:D+1),0);
psi2 = rho^2*k(x_star,x2,ModelInfo.hyp(1:D+1),0) + k(x_star,x2,ModelInfo.hyp(D+2:2*D+2),0);
psi = [psi1 psi2];

% calculate prediction
f = psi*(L'\(L\y));

v = rho^2*k(x_star, x_star, ModelInfo.hyp(1:D+1),0) ...
    + k(x_star, x_star, ModelInfo.hyp(D+2:2*D+2),0) ...
    - psi*(L'\(L\psi'));

v = diag(v);
