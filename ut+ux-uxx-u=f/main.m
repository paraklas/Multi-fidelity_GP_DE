clear all
clc
close all

addpath ../Utilities

rng('default')

global ModelInfo

%% Setup
Ntr = [30 10];
dim = 2;
Nbc = [3 4];
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

% Initial/Boundary conditions
lt = bsxfun(@plus,lb(1),bsxfun(@times,   lhsdesign(Nbc(1),1)    ,(ub(1)-lb(1))));
onet = ones(Nbc(1),1);
lx = bsxfun(@plus,lb(2),bsxfun(@times,   lhsdesign(Nbc(2),1)    ,(ub(2)-lb(2))));
onex = ones(Nbc(2),1);
b1 = [lt onet]; b2 = [lt 0.0*onet]; b3 = [0.0*onex lx];
ModelInfo.x0 = [b1; b2; b3];
ModelInfo.u0 = Exact_solution(ModelInfo.x0);
ModelInfo.u0 = ModelInfo.u0 + 0.02*randn(size(ModelInfo.u0));

ModelInfo.hyp = log([1 ones(1,dim) 1 ones(1,dim) exp(2) 10^-3 10^-3 10^-3]);
[ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -5000);

ModelInfo.hyp(end-3)
exp(ModelInfo.hyp)


[NegLnLike]=likelihood(ModelInfo.hyp);

 
%% Make Predictions & Plot results
save_plots = 1;
time_dep = 1;
PlotResults