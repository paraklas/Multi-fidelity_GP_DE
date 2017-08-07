clear all
clc
close all

addpath ../Utilities

rng('default')

global ModelInfo


%% Setup
Ntr = [0 3];
dim = 1;
Nbc = 1*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = eps;
ModelInfo.jitter=jitter;

%% Optimize model

% Training data for RHS
ModelInfo.x2 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(2),dim)    ,(ub-lb)));
ModelInfo.y2 = RHS_H(ModelInfo.x2);
ModelInfo.y2 = ModelInfo.y2 + 0.05*randn(size(ModelInfo.y2));

ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));
ModelInfo.y1 = RHS_L(ModelInfo.x1);
ModelInfo.y1 = ModelInfo.y1 + 0.3*randn(size(ModelInfo.y1));

ModelInfo.x0 = 0.1;
ModelInfo.u0 = Exact_solution(ModelInfo.x0);

% hyp = [sigma1 theta1 sigma2 theta2 rho sigma_n]
ModelInfo.hyp = log([1 1 1 1 exp(2) 10^-3 10^-6 0]);
[ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -5000);

ModelInfo.hyp(end-3)
exp([ModelInfo.hyp(1:end-4) ModelInfo.hyp(end-2:end)])
[NLML,~]=likelihood(ModelInfo.hyp);
fprintf(1,'Negative Log Marginal likelihood: %e\n', NLML);


%% Make Predictions & Plot results
save_plots = 1;
PlotResults