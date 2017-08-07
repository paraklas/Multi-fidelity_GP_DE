clear all
clc
close all

set(0,'defaulttextinterpreter','latex')
color1 = [0 90 178]/255;
color2 = [208 61 3]/255;

addpath ../Utilities

rng('default')

global ModelInfo

%% Setup
Ntr = 4;
dim = 2;
Nbc = 25*ones(1,dim);
lb = zeros(1,dim);
ub = ones(1,dim);
jitter = 1e-12;
ModelInfo.jitter = jitter;
noise = 0.0;

iter = 0;
max_iter = 80;
plt = 1;

%% Optimize model

% Training data for RHS
ModelInfo.x1 = bsxfun(@plus,lb,bsxfun(@times,   lhsdesign(Ntr(1),dim)    ,(ub-lb)));
%ModelInfo.x1 = [0.25 0.25; 0.25 0.75; 0.75 0.25; 0.75 0.75];
ModelInfo.y1 = RHS(ModelInfo.x1);
ModelInfo.y1 = ModelInfo.y1 + noise*randn(size(ModelInfo.y1));

% Boundary conditions
xbc = GetBoundary(dim, Nbc);
ModelInfo.x0 = bsxfun(@plus,lb,bsxfun(@times,xbc,(ub-lb)));
ModelInfo.u0 = Exact_solution(ModelInfo.x0);
ModelInfo.u0 = ModelInfo.u0 + noise*randn(size(ModelInfo.u0));


% Test points
nn = 40;
x1 = linspace(lb(1), ub(1), nn)';
x2 = linspace(lb(2), ub(2), nn)';
[Xplot, Yplot] = meshgrid(x1,x2);
x_star = reshape([Xplot Yplot], nn^2, 2);
Nts = size(x_star,1);
    
Exact = Exact_solution(x_star);
Exact_rhs = RHS(x_star);
Vx_star = zeros(Nts,Nts);

Exactplot = griddata(x_star(:,1),x_star(:,2),Exact,Xplot,Yplot,'cubic');
Exactplot_rhs = griddata(x_star(:,1),x_star(:,2),Exact_rhs,Xplot,Yplot,'cubic');

% Initialize hyp := [sigma1 theta1 sigma_n0 sigma_n1]
ModelInfo.hyp = log([1 ones(1,dim)]);
% [hyp_init, ~] = ga(@likelihood_init,dim+3,[],[],[],[],-3*ones(1,dim+3),3*ones(1,dim+3),[]);
% ModelInfo.hyp = hyp_init;

fig = figure('units','normalized','outerposition',[0 0 1 1]);
% movObj = QTWriter('active_learning.mov');
movObj = QTWriter('active_learning_fast.mov');

cap=imread('caption.png');

% error = zeros(max_iter,1);
% error_rhs = zeros(max_iter,1);
hyps = zeros(max_iter, 3);
while (iter <= max_iter)

    iter = iter + 1;
    clf
    
    [ModelInfo.hyp,~,~] = minimize(ModelInfo.hyp, @likelihood, -2000);       
    [NLML, ~]=likelihood(ModelInfo.hyp);
    fprintf(1,'Negative Log Marginal likelihood: %e\n', NLML);
    exp(ModelInfo.hyp)
    hyps(iter,:) = exp(ModelInfo.hyp);

    % Prediction
    [Kpred, Kvar] = predictor(x_star);
    [Kpred_rhs, Kvar_rhs] = predictor_rhs(x_star);

    % Active learning acquisition (IMSE)
%     for i = 1:Nts
%         Vx_star(i,:) = Kvar_rhs - cvar_rhs(x_star, x_star(i,:));
%     end
%     [~, i] = max(sum(Vx_star,1));
%     new_X = x_star(i,:);

%     Active learning acquisition (Variance)
%     [v1, i1] = max(Kvar);
%     [v2, i2] = max(Kvar_rhs);
%     if (v1/norm(Kvar)) > (v2/norm(Kvar_rhs))
%         new_X = x_star(i1,:);
%     else
%         new_X = x_star(i2,:);
%     end
    [~, i] = max(Kvar_rhs);
    new_X = x_star(i,:);
%    i = randi(Nts);
  

    error(iter) = norm(Kpred-Exact,2)/norm(Exact,2);
    error_rhs(iter) = norm(Kpred_rhs-Exact_rhs,2)/norm(Exact_rhs,2);

    fprintf(1,'Iteration: %d, error_u = %e, error_f = %e\n', iter, error(iter), error_rhs(iter));
    if plt == 1
        
        Predplot = griddata(x_star(:,1),x_star(:,2),Kpred,Xplot,Yplot,'cubic');
        Varplot = griddata(x_star(:,1),x_star(:,2),Kvar,Xplot,Yplot,'cubic');
        Predplot_rhs = griddata(x_star(:,1),x_star(:,2),Kpred_rhs,Xplot,Yplot,'cubic');
        Varplot_rhs = griddata(x_star(:,1),x_star(:,2),Kvar_rhs,Xplot,Yplot,'cubic');
    
        subplot(2,4,1);
        hold
        pcolor(Xplot,Yplot,real(sqrt(Varplot)));
        shading interp
        plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',10);
        %plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','m', 'MarkerEdgeColor','m','MarkerSize',10);
        plot(new_X(1), new_X(2),'p','MarkerFaceColor','k', 'MarkerEdgeColor','b','MarkerSize',18);
        colorbar
        title('$u(x)$ std');
        xlabel('$x_{1}$')
        ylabel('$x_{2}$')
        axis square
        set(gca,'FontSize',18);
        set(gcf, 'Color', 'w');
        
        subplot(2,4,2);
        hold
        pcolor(Xplot,Yplot,abs(Exactplot-Predplot)/max(abs(Exact)));
        shading interp
        plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',10);
        %plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','m', 'MarkerEdgeColor','m','MarkerSize',10);
        plot(new_X(1), new_X(2),'p','MarkerFaceColor','k', 'MarkerEdgeColor','b','MarkerSize',18);
        colorbar
        title('$u(x)$ error');
        xlabel('$x_{1}$')
        ylabel('$x_{2}$')
        axis square
        set(gca,'FontSize',18);
        set(gcf, 'Color', 'w');
        
        xxx = subplot(2,4,5);
        hold
        pcolor(Xplot,Yplot,real(sqrt(Varplot_rhs)));
        colormap(xxx,cool);
        shading interp
        plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',10);
        %plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','m', 'MarkerEdgeColor','m','MarkerSize',10);
        plot(new_X(1), new_X(2),'p','MarkerFaceColor','k', 'MarkerEdgeColor','b','MarkerSize',18);
        colorbar
        title('$f(x)$ std');
        xlabel('$x_{1}$')
        ylabel('$x_{2}$')
        axis square
        set(gca,'FontSize',18);
        set(gcf, 'Color', 'w');
        
        
        subplot(2,4,6);
        hold
        pcolor(Xplot,Yplot,abs(Exactplot_rhs-Predplot_rhs)/max(abs(Exact_rhs)));
        shading interp
        plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',10);
        %plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','m', 'MarkerEdgeColor','m','MarkerSize',10);
        plot(new_X(1), new_X(2),'p','MarkerFaceColor','k', 'MarkerEdgeColor','b','MarkerSize',18);
        colorbar
        title('$f(x)$ error');
        xlabel('$x_{1}$')
        ylabel('$x_{2}$')
        axis square
        set(gca,'FontSize',18);
        set(gcf, 'Color', 'w');
        
        

        subplot(2,4,[3 4 7 8]);
        hold
        h(1) = plot(1:1:iter, error, 'Color', color2, 'LineWidth',3);
        h(2) = plot(1:1:iter, error_rhs, '--', 'Color', color1, 'LineWidth',3);
        xlabel('$n$')
        ylabel('Relative error')
        hl = legend(h, 'Solution', 'Forcing');%,'Location','southoutside');
        set(hl,'Interpreter','latex')
        set(gca,'FontSize',18);
        set(gcf, 'Color', 'w');
        xlim([0 90]);
        ylim([1e-10 1]);
        axis square
        % set(gca,'YTick',1e-8:1e-2:1);

        % legend('eu_imse','ef_imse','eu','ef');
        buft = sprintf('Iteration: %d, Number of training points: %d', iter-1, size(ModelInfo.x1,1));
        title(buft);
        set(gca,'YScale','log');
        legend boxoff
        
        
        axes('position',[0.18,-0.18,0.5,0.5])
        imshow(cap)
        
%         frame = getframe(fig);
%         im = frame2im(frame);
%         [imind, cm] = rgb2ind(im, 256);
%         if iter == 1
%             imwrite(imind,cm,'movie.gif','gif','Loopcount',Inf);
%         else
%             imwrite(imind,cm,'movie.gif','gif','WriteMode','append','DelayTime',2);
%         end
        
        % Vary the frame-rate
        movObj.FrameRate = 1.4;

        % Write each frame to the file
        writeMovie(movObj,getframe(fig));



        
%         w = waitforbuttonpress;
%         if w == 0
%             buf = sprintf('it_%d.png', iter);
%             export_fig(buf, '-r300');
%         end
        
        drawnow;
        
    end

    ModelInfo.x1 = [ModelInfo.x1; new_X];
    ModelInfo.y1 = [ModelInfo.y1; RHS(new_X) + noise*randn()];
    %ModelInfo.hyp = log([1 0.5*ones(1,dim) 10^-6 10^-6]);
%     ModelInfo.hyp = hyp_init;

end

% Set palindromic looping flag
movObj.Loop = 'backandforth';

% Finish writing movie and close file
close(movObj);

figure(2)
clf
hold
semilogy(error,'r')
semilogy(error_rhs,'b')
