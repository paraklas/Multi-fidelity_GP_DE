clear all
clc
close all

addpath ../Utilities/

rng('default')

global ModelInfo
ModelInfo.alpha = 0.3;

%% Setup
Ntr = [15 4];
N0 = 2;
dim = 1;
Nbc = 1*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = 1e-8;
ModelInfo.jitter = jitter;
        
%% Optimize model

% Training data for RHS
ModelInfo.x2 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(2),dim)    ,(ub-lb)));
ModelInfo.x0 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(N0,dim)    ,(ub-lb)));
ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));


ModelInfo.y1 = RHS_L(ModelInfo.x1);
ModelInfo.y2 = RHS_H(ModelInfo.x2);
% Add Noise 
ModelInfo.y1 = ModelInfo.y1 + 0.3*randn(size(ModelInfo.y1));
ModelInfo.y2 = ModelInfo.y2 + 0.05*randn(size(ModelInfo.y2));

ModelInfo.u0 = Exact_solution(ModelInfo.x0);

ModelInfo.hyp = log([1 ones(1,dim) 1 ones(1,dim) exp(1.5) 10^-3 10^-3 10^-3]);
[ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -5000);

ModelInfo.hyp(end-3)
exp(ModelInfo.hyp)


[NegLnLike]=likelihood(ModelInfo.hyp);

 
%% Make Predictions & Plot results
save_plots = 1;
PlotResults