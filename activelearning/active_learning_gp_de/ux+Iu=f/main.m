clear all
clc
close all

addpath ../Utilities

rng('default')

set(0,'defaulttextinterpreter','latex')

global ModelInfo


%% Setup
Ntr = 8;
Nts = 100;
dim = 1;
Nbc = 1*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = 1e-12;
noise = 0.0;
ModelInfo.jitter=jitter;

iter = 0;
max_iter = 20;
plt = 1;

%% Optimize model

% Training data for RHS
ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));
ModelInfo.y1 = RHS(ModelInfo.x1);
ModelInfo.y1 = ModelInfo.y1 + noise*randn(size(ModelInfo.y1));

% Traing data for u (bounday/initial data, etc.)
ModelInfo.x0 = 0.0;
ModelInfo.u0 = Exact_solution(ModelInfo.x0);

% Test points
x_star = linspace(lb,ub,Nts)';
Exact = Exact_solution(x_star);
Exact_rhs = RHS(x_star);
Vx_star = zeros(Nts,Nts);

% iInitialize hyp := [sigma1 theta1 sigma_n0 sigma_n1]
ModelInfo.hyp = log([0.2 0.2 10^-6 10^-6]);

figure('units','normalized','outerposition',[0 0 1 1])
color2 = [217,95,2]/255;
while (iter <= max_iter)

    iter = iter + 1;
    clf
    
    [ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -2000);
    [NLML,~]=likelihood(ModelInfo.hyp);
    fprintf(1,'Negative Log Marginal likelihood: %e\n', NLML);
    exp(ModelInfo.hyp)

    % Prediction
    [Kpred, Kvar] = predictor(x_star);
    [Kpred_rhs, Kvar_rhs] = predictor_rhs(x_star);

    % Active learning acquisition (IMSE)
    for i = 1:Nts
        Vx_star(i,:) = Kvar_rhs - cvar_rhs(x_star, x_star(i,:));
    end
    [~, i] = max(sum(Vx_star,1));
    new_X = x_star(i);

    % Active learning acquisition (Variance)
    % [~, i] = max(Kvar_rhs);
    % new_X = x_star(i);

    error = norm(Kpred-Exact,2)/norm(Exact,2);
    error_rhs = norm(Kpred_rhs-Exact_rhs,2)/norm(Exact_rhs,2);

    fprintf(1,'Iteration: %d, error_u = %e, error_f = %e, new_X = %f\n', iter, error, error_rhs, new_X);
    if plt == 1
        subplot(1,2,1);
        hold
        plot(x_star,Exact,'b','LineWidth',3);
        plot(x_star, Kpred,'r--','LineWidth',3);
        [l,p] = boundedline(x_star, Kpred, 2.0*sqrt(Kvar), ':', 'alpha','cmap', color2);
        outlinebounds(l,p);    
        yLimits = get(gca,'YLim');
        plot(new_X*ones(10,1), linspace(yLimits(1),yLimits(2),10),'k--','LineWidth',2);
        xlabel('$x$')
        ylabel('$u(x)$')
        axis square
        ylim([yLimits(1) yLimits(2)]);
        set(gca,'FontSize',24);
        set(gcf, 'Color', 'w');
        
        
        subplot(1,2,2);
        hold
        plot(x_star, Exact_rhs,'b','LineWidth',3);
        plot(x_star,Kpred_rhs,'r--','LineWidth',3);
        [l,p] = boundedline(x_star, Kpred_rhs, 2.0*sqrt(Kvar_rhs), ':', 'alpha','cmap', color2);
        outlinebounds(l,p);            
        plot(ModelInfo.x1, ModelInfo.y1,'kx','MarkerSize',16, 'LineWidth',3);
        yLimits = get(gca,'YLim');
        plot(new_X*ones(10,1), linspace(yLimits(1),yLimits(2),10),'k--','LineWidth',2);        xlabel('$x$')
        ylabel('$f(x)$')
        axis square
        ylim([yLimits(1) yLimits(2)]);
        set(gca,'FontSize',24);
        set(gcf, 'Color', 'w');
    
        w = waitforbuttonpress;
        
    end

    ModelInfo.x1 = [ModelInfo.x1; new_X];
    ModelInfo.y1 = [ModelInfo.y1; RHS(new_X) + noise*randn()];
    hyp = ModelInfo.hyp;

end
