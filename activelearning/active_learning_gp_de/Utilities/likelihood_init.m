function [NLML]=likelihood_init(hyp)

global ModelInfo
x0 = ModelInfo.x0;
x1 = ModelInfo.x1;

u0 = ModelInfo.u0;
y1 = ModelInfo.y1;
y=[u0;y1];

jitter = ModelInfo.jitter;

sigma_n1 = exp(hyp(end));
sigma_n0 = exp(hyp(end-1));

n0 = size(x0,1);
[n1,D] = size(x1);
n = n0+n1;

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

ModelInfo.L = L;

if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

alpha = L'\(L\y);
NLML = 0.5*y'*alpha + sum(log(diag(L))) + log(2*pi)*n/2;


% D_NLML = 0*hyp;
% Q =  L'\(L\eye(n)) - alpha*alpha';
% for i=1:D+1
%     DK00 = g(x0, x0, hyp(1:D+1),i);
%     DK01 = h(x0, x1, hyp(1:D+1),i);
%     DK10 = DK01';
%     DK11 = k(x1, x1, hyp(1:D+1),i);
%     
%     DK = [DK00 DK01;
%         DK10 DK11];
%     D_NLML(i) = sum(sum(Q.*DK))/2;
% end
% 
% D_NLML(end-1) = sigma_n0*trace(Q(1:n0,1:n0))/2;
% D_NLML(end) = sigma_n1*trace(Q(n0+1:n0+n1,n0+1:n0+n1))/2;