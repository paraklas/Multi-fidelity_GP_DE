clear all
clc
close all

addpath ../Utilities

rng('default')

global ModelInfo

%% Setup
Ntr = [40 10];
dim = 2;
Nbc = 5*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = eps;
ModelInfo.jitter = jitter;

%% Optimize model

% Training data for RHS

ModelInfo.x2 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(2),dim)    ,(ub-lb)));
ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));

ModelInfo.y1 = RHS_L(ModelInfo.x1);
ModelInfo.y2 = RHS_H(ModelInfo.x2);
 
ModelInfo.y1 = ModelInfo.y1 + 0.3*randn(size(ModelInfo.y1));
ModelInfo.y2 = ModelInfo.y2 + 0.05*randn(size(ModelInfo.y2));

% Boundary conditions
xbc = GetBoundary(dim, Nbc);


ModelInfo.x0 = bsxfun(@plus,lb,bsxfun(@times,xbc,(ub-lb)));
ModelInfo.u0 = Exact_solution(ModelInfo.x0);
ModelInfo.u0 = ModelInfo.u0 + 0.01*randn(size(ModelInfo.u0));

ModelInfo.hyp = log([1 ones(1,dim) 1 ones(1,dim) exp(1) 10^-3 10^-3 10^-3]);
tic
[ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -2000);
toc

ModelInfo.hyp(end-3)
exp(ModelInfo.hyp)

[NegLnLike]=likelihood(ModelInfo.hyp);


%% Make Predictions & Plot results
save_plots = 1;
time_dep = 0;
PlotResults