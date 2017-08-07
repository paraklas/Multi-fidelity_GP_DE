function [v] = cvar_rhs(x_ss, x_s)

global ModelInfo

hyp = ModelInfo.hyp;
x0 = ModelInfo.x0;
x1 = [ModelInfo.x1; x_s];

jitter = ModelInfo.jitter;

sigma_n0 = exp(hyp(end-1));
sigma_n1 = exp(hyp(end));

[n1,~] = size(x1);
[n0,D] = size(x0);
n = n0 + n1;

K00 = g(x0, x0, hyp(1:D+1),0);

K01 = h(x0, x1, hyp(1:D+1),0);

K10 = K01';

K11 = k(x1, x1, hyp(1:D+1),0);

K00 = K00 + eye(n0).*sigma_n0;
K11 = K11 + eye(n1).*sigma_n1;

K = [K00 K01;
    K10 K11];

K = K + eye(n).*jitter;

% Cholesky factorisation
[L,p]=chol(K,'lower');
if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

psi1 = h(x0, x_ss, hyp(1:D+1),0)';
psi2 = k(x1, x_ss, hyp(1:D+1),0)';
psi = [psi1 psi2];

% calculate prediction
v = k(x_ss, x_ss, ModelInfo.hyp(1:D+1),0) ...
    - psi*(L'\(L\psi'));

v = diag(v);

end