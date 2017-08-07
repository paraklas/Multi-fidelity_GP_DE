function [NLML,D_NLML]=likelihood(hyp)

global ModelInfo
x0 = ModelInfo.x0;
x1 = ModelInfo.x1;
x2 = ModelInfo.x2;

u0 = ModelInfo.u0;
y1 = ModelInfo.y1;
y2 = ModelInfo.y2;
y=[u0;y1;y2];

jitter = ModelInfo.jitter;

sigma_n2 = exp(hyp(end));
sigma_n1 = exp(hyp(end-1));
sigma_n0 = exp(hyp(end-2));
rho = hyp(end-3);

n0 = size(x0,1);
n1 = size(x1,1);
[n2,D] = size(x2);
n = n0+n1+n2;

K00 = rho^2*g(x0, x0, hyp(1:D+1),0) ...
    + g(x0, x0, hyp(D+2:2*D+2),0);

K01 = rho*h(x0, x1, hyp(1:D+1),0);
K02 = rho^2*h(x0, x2, hyp(1:D+1),0) ...
    + h(x0, x2, hyp(D+2:2*D+2),0);

K10 = K01';
K20 = K02';

K11 = k(x1, x1, hyp(1:D+1),0);
K12 = rho*k(x1, x2, hyp(1:D+1),0);
K21 = rho*k(x2, x1, hyp(1:D+1),0);
K22 = (rho^2)*k(x2, x2, hyp(1:D+1),0) + k(x2, x2, hyp(D+2:2*D+2),0);

K00 = K00 + eye(n0).*sigma_n0;
K11 = K11 + eye(n1).*sigma_n1;
K22 = K22 + eye(n2).*sigma_n2;

K = [K00 K01 K02;
    K10 K11 K12;
    K20 K21 K22];

K = K + eye(n).*jitter;

% Cholesky factorisation
[L,p]=chol(K,'lower');

ModelInfo.L = L;

if p > 0
    fprintf(1,'Covariance is ill-conditioned\n');
end

alpha = L'\(L\y);
NLML = 0.5*y'*alpha + sum(log(diag(L))) + log(2*pi)*n/2;


D_NLML = 0*hyp;
Q =  L'\(L\eye(n)) - alpha*alpha';
for i=1:D+1
    DK00 = rho^2*g(x0, x0, hyp(1:D+1),i);
    DK01 = rho*h(x0, x1, hyp(1:D+1),i);
    DK02 = rho^2*h(x0, x2, hyp(1:D+1),i);
    DK10 = DK01';
    DK20 = DK02';
    
    DK11 = k(x1, x1, hyp(1:D+1),i);
    DK12 = rho*k(x1, x2, hyp(1:D+1),i);
    DK21 = rho*k(x2, x1, hyp(1:D+1),i);
    DK22 = rho^2*k(x2, x2, hyp(1:D+1),i);
    
    DK = [DK00 DK01 DK02;
        DK10 DK11 DK12;
        DK20 DK21 DK22];
    D_NLML(i) = sum(sum(Q.*DK))/2;
end

for i=D+2:2*D+2
    DK00 = g(x0, x0, hyp(D+2:2*D+2),i-(D+1));
    DK01 = zeros(n0,n1);
    DK02 = h(x0, x2, hyp(D+2:2*D+2),i-(D+1));
    DK10 = DK01';
    DK20 = DK02';
    
    DK11 = zeros(n1,n1);
    DK12 = zeros(n1,n2);
    DK21 = zeros(n2,n1);
    DK22 = k(x2, x2, hyp(D+2:2*D+2),i-(D+1));
    
    DK = [DK00 DK01 DK02;
        DK10 DK11 DK12;
        DK20 DK21 DK22];
    D_NLML(i) = sum(sum(Q.*DK))/2;
end

DK00 = 2*rho*g(x0, x0, hyp(1:D+1),0);
DK01 = h(x0, x1, hyp(1:D+1),0);
DK02 = 2*rho*h(x0, x2, hyp(1:D+1),0);
DK10 = DK01';
DK20 = DK02';

DK11 = zeros(n1,n1);
DK12 = k(x1, x2, hyp(1:D+1),0);
DK21 = k(x2, x1, hyp(1:D+1),0);
DK22 = (2*rho)*k(x2, x2, hyp(1:D+1),0);
DK = [DK00 DK01 DK02;
    DK10 DK11 DK12;
    DK20 DK21 DK22];

D_NLML(end-3) = sum(sum(Q.*DK))/2;

D_NLML(end-2) = sigma_n0*trace(Q(1:n0,1:n0))/2;
D_NLML(end-1) = sigma_n1*trace(Q(n0+1:n0+n1,n0+1:n0+n1))/2;
D_NLML(end) = sigma_n2*trace(Q(n0+n1+1:end,n0+n1+1:end))/2;