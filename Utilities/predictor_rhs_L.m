function [f, v] = predictor_rhs_L(x_star)

global ModelInfo
hyp = ModelInfo.hyp;
x1 = ModelInfo.x1;

y1 = ModelInfo.y1;
y=y1;

jitter = ModelInfo.jitter;

sigma_n1 = exp(hyp(end-1));

[n1,D] = size(x1);
n = n1;

K = k(x1, x1, hyp(1:D+1),0);

K = K + eye(n1).*sigma_n1;

K = K + eye(n).*jitter;

% Cholesky factorisation
[L,p]=chol(K,'lower');

if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

psi = k(x_star,x1,ModelInfo.hyp(1:D+1),0);

% calculate prediction
f = psi*(L'\(L\y));

v = k(x_star, x_star, ModelInfo.hyp(1:D+1),0) - psi*(L'\(L\psi'));

v = diag(v);
