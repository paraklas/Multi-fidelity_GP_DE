clear all
clc
close all

rng('default')

global ModelInfo


%% Setup
dim = 1;
Nbc = 1*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = eps;
ModelInfo.jitter=jitter;

nhigh = 3:1:7;
nlow = 0:1:30;


[X, Y] = meshgrid(nlow,nhigh);
pp = reshape([X Y], length(nlow)*length(nhigh), 2);

error_u = zeros(length(pp),1);
error_rhs = zeros(length(pp),1);

for i = 1:length(pp)
    fprintf(1,'Iteration %d/%d', i, length(pp));
    Ntr = pp(i,:);
    
    %% Optimize model
    rng('default')
    % Training data for RHS
    %     ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));
    ModelInfo.x2 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(2),dim)    ,(ub-lb)));
    ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));
    
    %     ModelInfo.x1 = linspace(0,1,Ntr(1))';
    %     ModelInfo.x2 = linspace(0,1,Ntr(2))';
    ModelInfo.y1 = RHS_L(ModelInfo.x1);
    ModelInfo.y2 = RHS_H(ModelInfo.x2);
    
    ModelInfo.y1 = ModelInfo.y1 + 0.001*randn(size(ModelInfo.y1));
    ModelInfo.y2 = ModelInfo.y2 + 0.00001*randn(size(ModelInfo.y2));
    
    ModelInfo.x0 = 0.0;
    ModelInfo.u0 = Exact_solution(ModelInfo.x0);
    
    % hyp = [sigma1 theta1 sigma2 theta2 rho sigma_n]
    ModelInfo.hyp = log([1 1 1 1 exp(1.3) 10^-3 10^-3 10^-3]);
    [ModelInfo.hyp,~,~] = minimize_nv(ModelInfo.hyp, @likelihood, -200);
    % ModelInfo.hyp(end-1) = 0.0;
    %ModelInfo.hyp = fminunc(@likelihood,ModelInfo.hyp);
    ModelInfo.hyp(end-3)
    [NLML,~]=likelihood(ModelInfo.hyp);
    fprintf(1,'Negative Log Marginal likelihood: %e\n', NLML);
    
    %% Make Predictions & Plot results
    nn = 200;
    Xts = linspace(lb(1), ub(1), nn)';
    Nts = size(Xts,1);
    
    [Kpred, Kvar] = predictor(Xts);
    [Kpred_rhs_H, Kvar_rhs_H] = predictor_rhs_H(Xts);
    
    Exact = Exact_solution(Xts);
    Exact_rhs = RHS_H(Xts);

    error_u(i) = (norm(Kpred-Exact,2)/norm(Exact,2));
    disp(error_u(i))
    error_rhs(i) = (norm(Kpred_rhs_H-Exact_rhs,2)/norm(Exact_rhs,2));
    disp(error_rhs(i))
end