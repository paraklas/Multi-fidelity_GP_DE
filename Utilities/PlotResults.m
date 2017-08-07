set(0,'defaulttextinterpreter','latex')

if dim == 1
    nn = 200;
    Xts = linspace(lb(1), ub(1), nn)';
    Nts = size(Xts,1);
elseif dim == 2
    nn = 50;
    x1 = linspace(lb(1), ub(1), nn)';
    x2 = linspace(lb(2), ub(2), nn)';
    [Xplot, Yplot] = meshgrid(x1,x2);
    Xts = reshape([Xplot Yplot], nn^2, 2);
    Nts = size(Xts,1);
else
    rng('default');
    Nts = 1000;
    xts = rand(Nts,dim);
    Xts = bsxfun(@plus,lb,bsxfun(@times,xts,(ub-lb)));
    idx = 1;
    nplt = 80;
    Xplot = [linspace(lb(idx), ub(idx),nplt)' ones(nplt,dim-1)/4];
end

% Predictor
[Kpred, Kvar] = predictor(Xts);
[Kpred_rhs_H, Kvar_rhs_H] = predictor_rhs_H(Xts);
[Kpred_rhs_L, Kvar_rhs_L] = predictor_rhs_L(Xts);

Exact = Exact_solution(Xts);
Exact_rhs_H = RHS_H(Xts);
Exact_rhs_L = RHS_L(Xts);

fprintf(1,'Relative L2 error u: %e\n', (norm(Kpred-Exact,2)/norm(Exact,2)));
fprintf(1,'Relative L2 error RHS_H: %e\n', (norm(Kpred_rhs_H-Exact_rhs_H,2)/norm(Exact_rhs_H,2)));
fprintf(1,'Relative L2 error RHS_L: %e\n', (norm(Kpred_rhs_L-Exact_rhs_L,2)/norm(Exact_rhs_L,2)));


% buf = sprintf('Exact_d%d_n1-%d_n2-%d.txt', size(Xts,2), size(ModelInfo.x1,1), size(ModelInfo.x2,1));
% save(buf,'-ascii','-double','Exact');
% buf = sprintf('Kpred_d%d_n1-%d_n2-%d.txt', size(Xts,2), size(ModelInfo.x1,1), size(ModelInfo.x2,1));
% save(buf,'-ascii','-double','Kpred');


color2 = [217,95,2]/255;

if dim == 1
    fig = figure(1);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    %figure('name','1', 'units','normalized','outerposition',[0 0 1 1])
    clf
    subplot(1,2,1)
    clear h
    clear leg
    hold
    h(1) = plot(Xts, Exact,'k','LineWidth',3);
    h(2) = plot(Xts,Kpred,'r--','LineWidth',3);
    [l,h(3)] = boundedline(Xts, Kpred, 2.0*sqrt(Kvar), ':', 'alpha','cmap', color2);
    outlinebounds(l,h(3));
    h(4) = plot(ModelInfo.x0, ModelInfo.u0,'s', 'MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',12);
    hl = legend(h,'Exact solution', 'GP posterior mean', 'Two standard deviations band','Anchor point(s)','Location','southoutside');
    legend boxoff
    set(hl,'Interpreter','latex')
    xlabel('$x$')
    ylabel('$u(x)$')
%     if Ntr(1) > 0
%         title('(a)')
%     else
%         title('(c)')
%     end
    axis square
    set(gca,'FontSize',32);
    set(gcf, 'Color', 'w');
    
    
    subplot(1,2,2)
    clear h
    clear leg
    hold
    if Ntr(1) > 0
        h(1) = plot(Xts, Exact_rhs_H,'k','LineWidth',3);
        %h(2) = plot(Xts,Kpred_rhs_H,'r--','LineWidth',3);
        %[l,h(3)] = boundedline(Xts, Kpred_rhs_H, 2.0*sqrt(Kvar_rhs_H), ':', 'alpha','cmap', color2);
        %outlinebounds(l,h(3));
        
        h(2) = plot(Xts, Exact_rhs_L,'b--','LineWidth',1);
        %h(5) = plot(Xts,Kpred_rhs_L,'m--','LineWidth',2);
        %[l,h(6)] = boundedline(Xts, Kpred_rhs_L, 2.0*sqrt(Kvar_rhs_L), ':', 'alpha','cmap', color2);
        %outlinebounds(l,h(6));
        
        h(3) = plot(ModelInfo.x2, ModelInfo.y2,'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',12);
        
        
        h(4) = plot(ModelInfo.x1, ModelInfo.y1,'bo','MarkerSize',8);
    else
        h(1) = plot(Xts, Exact_rhs_H,'k','LineWidth',3);
        %h(2) = plot(Xts,Kpred_rhs_H,'r--','LineWidth',3);
        %[l,h(3)] = boundedline(Xts, Kpred_rhs_H, 2.0*sqrt(Kvar_rhs_H), ':', 'alpha','cmap', color2);
        %outlinebounds(l,h(3));
        
        %    h(2) = plot(Xts, Exact_rhs_L,'b--','LineWidth',3);
        %h(5) = plot(Xts,Kpred_rhs_L,'m--','LineWidth',2);
        %[l,h(6)] = boundedline(Xts, Kpred_rhs_L, 2.0*sqrt(Kvar_rhs_L), ':', 'alpha','cmap', color2);
        %outlinebounds(l,h(6));
        
        h(2) = plot(ModelInfo.x2, ModelInfo.y2,'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',12);
        
        
        %        h(4) = plot(ModelInfo.x1, ModelInfo.y1,'o', 'MarkerFaceColor','b', 'MarkerEdgeColor','b','MarkerSize',8);
        
    end
    
    
    buf1 = sprintf('High-fidelity training data (%d points)', Ntr(2));
    buf2 = sprintf('Low-fidelity training data (%d points)', Ntr(1));
    if Ntr(1) > 0
        leg{1} = 'Exact hight-fidelity forcing';% leg{2} ='Predictive mean (high-fidelity)';
        %leg{3} ='Two standard deviation band (high-fidelity)';
        
        leg{2} = 'Exact low-fidelity forcing';% leg{5} ='Predictive mean (low-fidelity)';
        %leg{6} ='Two standard deviation band (low-fidelity)';
        
        leg{3} = buf1;
        
        leg{4} = buf2;
    else
        leg{1} = 'Exact hight-fidelity forcing';% leg{2} ='Predictive mean (high-fidelity)';
        %leg{3} ='Two standard deviation band (high-fidelity)';
        
        %leg{2} = 'Exact low-fidelity forcing';% leg{5} ='Predictive mean (low-fidelity)';
        %leg{6} ='Two standard deviation band (low-fidelity)';
        
        leg{2} = buf1;
        
        %    leg{4} = buf2;
        
    end
    
    xlabel('$x$')
    if Ntr(1) > 0
        ylabel('$f_{1}(x), f_{2}(x)$')
    else
        ylabel('$f(x)$')
    end
    hl = legend(h,leg,'Location','southoutside'); legend boxoff
    set(hl,'Interpreter','latex')
%     if Ntr(1) > 0
%         title('(b)')
%     else
%         title('(d)')
%     end

    axis square
    set(gca,'FontSize',32);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig solution.png -r300
    end
    
    if Ntr(1) > 0
        fig = figure(2);
        set(fig,'units','normalized','outerposition',[0 0 1 1])
        clf
        hold
        h(1) = plot(Exact_rhs_L, Exact_rhs_H,'b','LineWidth',3);
        xlabel('$f_{1}(x)$');
        ylabel('$f_{2}(x)$');
        title('$f_{1}(x),f_{2}(x)$ Cross-correlation');
        axis square
        set(gca,'FontSize',48);
        set(gca,'Xtick',[]);
        set(gca,'Ytick',[]);
        set(gcf, 'Color', 'w');
        if save_plots  == 1
            export_fig correlation.png -r300
        end
    end
    
elseif dim > 2
    miny = min(Exact);
    maxy = max(Exact);
    
    [ev, idx] = sort(Exact);
    
    fig = figure(1);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    hold
    plot(ev, Kpred(idx), 'ko','MarkerSize', 10);
    plot(ev, ev, 'r-.','LineWidth',4);
    
    xlim([miny maxy]);
    ylim([miny maxy]);
    axis square
    buf = {'GP posterior mean','Exact'};
    
    xlabel('$u(x)$')
    ylabel('$\overline{u}(x)$')
    hl = legend(buf,'Location','southoutside'); legend boxoff
    set(hl,'Interpreter','latex')
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig scatter.png -r300
    end
    
    % 1D Predictor
    Eplt = Exact_solution(Xplot);
    [Kp_plt, Kv_plt] = predictor(Xplot);
    
    
    fig = figure(2);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    hold
    plot(Xplot(:,1), Eplt,'k','LineWidth',3)
    plot(Xplot(:,1),Kp_plt,'r--','LineWidth',3)
    [l,p] = boundedline(Xplot(:,1), Kp_plt, 2.0*sqrt(Kv_plt), ':', 'alpha','cmap', color2);
    outlinebounds(l,p);
    hl = legend('Exact solution', 'GP posterior mean', 'Two standard deviation band','Location','eastoutside');
    legend boxoff
    set(hl,'Interpreter','latex')
    xlabel('$x_1$')
    ylabel('$u(x)$')
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig solution.png -r300
    end
    
    Rhs_H_plt = RHS_H(Xplot);
    Rhs_L_plt = RHS_L(Xplot);
    [Rhs_H_p_plt, Rhs_H_v_plt] = predictor_rhs_H(Xplot);
    [Rhs_L_p_plt, Rhs_L_v_plt] = predictor_rhs_L(Xplot);
    
    fig = figure(3);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    hold
    h(1) = plot(Xplot(:,1), Rhs_H_plt,'b','LineWidth',3);
    h(2) = plot(Xplot(:,1),Rhs_H_p_plt,'r--','LineWidth',3);
    [l,h(3)] = boundedline(Xplot(:,1), Rhs_H_p_plt, 2.0*sqrt(Rhs_H_v_plt), ':', 'alpha','cmap', color2);
    outlinebounds(l,h(3));
    
    h(4) = plot(Xplot(:,1), Rhs_L_plt,'k','LineWidth',2);
    h(5) = plot(Xplot(:,1),Rhs_L_p_plt,'m--','LineWidth',2);
    [l,h(6)] = boundedline(Xplot(:,1), Rhs_L_p_plt, 2.0*sqrt(Rhs_L_v_plt), ':', 'alpha','cmap', color2);
    outlinebounds(l,h(6));
    
    
    buf1 = sprintf('High-fidelity training data (%d points)', Ntr(2));
    buf2 = sprintf('Low-fidelity training data (%d points)', Ntr(1));
    leg{1} = 'Exact hight-fidelity forcing'; leg{2} ='Predictive mean (high-fidelity)';
    leg{3} ='Two standard deviation band (high-fidelity)';
    
    leg{4} = 'Exact low-fidelity forcing'; leg{5} ='Predictive mean (low-fidelity)';
    leg{6} ='Two standard deviation band (low-fidelity)';
    
    xlabel('$x_1$')
    ylabel('$f_{1}(x), f_{2}(x)$')
    hl = legend(h,leg,'Location','eastoutside'); legend boxoff
    set(hl,'Interpreter','latex')
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig forcing.png -r300
    end
    
    
    fig = figure(4);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    clear leg
    clear h
    hold all
    h(1) = histogram(Exact,30,'Normalization','probability');
    h(2) = histogram(Kpred,30,'Normalization','probability');
    leg{1} = 'Exact solution $u(x)$';
    leg{2} = 'Predictive mean $\overline{u}(x)$';
    
    xlabel('$u(x), \overline{u}(x)$')
    ylabel('Normalized frequency')
    hl = legend(h,leg,'Location','southoutside'); legend boxoff
    set(hl,'Interpreter','latex')
    axis square
    xlim([-1.2,1.2])
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig histogram.png -r300
    end
    
    D = size(ModelInfo.x1,2);
    fig = figure(5);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    hold
    temp = [1./exp(ModelInfo.hyp(2:D+1));1./exp(ModelInfo.hyp(D+3:2*D+2))]';
    h = bar(temp);
    set(gca,'XTick',1:D)
    xlim([0,D+1])
    xlabel('Dimension $d=1,\ldots,10$')
    ylabel('ARD weights')
    hl = legend(h,{'High-fidelity','Low-fidelity'},'Location','southoutside');
    legend boxoff
    set(hl,'Interpreter','latex')
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig ARDweights.png -r300
    end
    
    
elseif dim == 2
    
    Exactplot = griddata(Xts(:,1),Xts(:,2),Exact,Xplot,Yplot,'cubic');
    Predplot = griddata(Xts(:,1),Xts(:,2),Kpred,Xplot,Yplot,'cubic');
    Varplot = griddata(Xts(:,1),Xts(:,2),Kvar,Xplot,Yplot,'cubic');
    
    Exactplot_rhs = griddata(Xts(:,1),Xts(:,2),Exact_rhs_H,Xplot,Yplot,'cubic');
    Predplot_rhs = griddata(Xts(:,1),Xts(:,2),Kpred_rhs_H,Xplot,Yplot,'cubic');
    Varplot_rhs = griddata(Xts(:,1),Xts(:,2),Kvar_rhs_H,Xplot,Yplot,'cubic');
    
    EH = griddata(Xts(:,1),Xts(:,2),Exact_rhs_H,Xplot,Yplot,'cubic');
    EL = griddata(Xts(:,1),Xts(:,2),Exact_rhs_L,Xplot,Yplot,'cubic');
    
    
    fig = figure(1);
    clear h;
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clf
    subplot(1,2,1);
    hold
    
    s1 = surf(Xplot,Yplot,Predplot);
    set(s1,'FaceColor',[100,100,100]/255,'EdgeColor','none', 'LineWidth',0.001,'FaceAlpha',0.35);
    material dull
    camlight
    lighting gouraud
    
    h(1) = plot3(ModelInfo.x0(:,1), ModelInfo.x0(:,2), ModelInfo.u0,'s', 'MarkerFaceColor','r', 'MarkerEdgeColor','r', 'MarkerSize',12);
    view(3)
    %zz = get(gca,'ZTick');
    %set(gca,'ZTickLabel',sprintf('%3.1f\n',zz));
    
    buf = cell(1,1);
    if time_dep == 1
        buf{1} = sprintf('Anchor points (%d points)',size(ModelInfo.x0,1));
    else
        buf{1} = sprintf('Anchor points (%d points)',size(ModelInfo.x0,1));
    end
    
    hl = legend(h,buf,'Location','northoutside');
    legend boxoff
    set(hl,'Interpreter','latex')
    
    if time_dep ==1
        xlabel('$t$')
        ylabel('$x$')
        zlabel('$\overline{u}(t,x)$')
    else
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$\overline{u}(t,x)$')
    end
    axis tight
    axis square
    set(gca,'FontSize',32);
    set(gcf, 'Color', 'w');
    
    clear h;
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,2);
    hold
    
    s1 = surf(Xplot,Yplot,Exactplot_rhs);
    set(s1,'FaceColor',[100,100,100]/255,'EdgeColor','none', 'LineWidth',0.001,'FaceAlpha',0.35);
    material dull    
    s2 = surf(Xplot,Yplot,EL);
    set(s2,'FaceColor','b','EdgeColor','none', 'LineWidth',0.001,'FaceAlpha',0.1);
    material dull
    camlight
    lighting gouraud
    
    h(1) = plot3(ModelInfo.x1(:,1), ModelInfo.x1(:,2), ModelInfo.y1,'bo','MarkerSize',8,'Linewidth',1);
    %     plot3(ModelInfo.x2(:,1), ModelInfo.x2(:,2), ModelInfo.y2+shift,'o', 'MarkerFaceColor',[30,30,30]/255, 'MarkerEdgeColor',[30,30,30]/255,'MarkerSize',10);
    h(2) = plot3(ModelInfo.x2(:,1), ModelInfo.x2(:,2), ModelInfo.y2,'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',12);
    
    view(3)
    %zz = get(gca,'ZTick');
    %set(gca,'ZTickLabel',sprintf('%3.1f\n',zz));
    
    buf{1} = sprintf('Low-fidelity training data (%d points)',size(ModelInfo.x1,1));
    buf{2} = sprintf('High-fidelity training data (%d points)',size(ModelInfo.x2,1));
    
    hl = legend(h,buf,'Location','northoutside'); legend boxoff
    set(hl,'Interpreter','latex')
    
    if time_dep == 1
        xlabel('$t$')
        ylabel('$x$')
        zlabel('$f_{1}(t,x), f_{2}(t,x)$')
    else
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$f_{1}(x), f_{2}(x)$')
    end
    axis tight
    axis square
    set(gca,'FontSize',32);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig solution.png -r300
    end
    
    
    
    fig = figure(2);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clear h
    clear buf
    clf
    hold
    contourf(Xplot,Yplot,sqrt(Varplot),16);
    shading interp
%     h(1) = plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',6);
%     h(2) = plot(ModelInfo.x2(:,1), ModelInfo.x2(:,2), 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',10);
%     h(3) = plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',10);
    
%     buf{1} = sprintf('Low-fidelity training data (%d points)',size(ModelInfo.x1,1));
%     buf{2} = sprintf('High-fidelity training data (%d points)',size(ModelInfo.x2,1));
%     buf{3} = sprintf('Anchor points (%d points)',size(ModelInfo.x0,1));
%     
%     hl = legend(h,buf,'Location','southoutside'); legend boxoff
%     set(hl,'Interpreter','latex')
    
%     title('Standard deviation');
    colorbar
    if time_dep == 1
        xlabel('$t$')
        ylabel('$x$')
    else
        xlabel('$x_1$')
        ylabel('$x_2$')
    end
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig variance.png -r300
    end
    
    fig = figure(3);
    set(fig,'units','normalized','outerposition',[0 0 1 1])
    clear h
    clear buf
    hold
    contourf(Xplot,Yplot,abs(Exactplot - Predplot)/max(max(abs(Exactplot))),14);
    shading interp
%     h(1) = plot(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',6);
%     h(2) = plot(ModelInfo.x2(:,1), ModelInfo.x2(:,2), 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',10);
%     h(3) = plot(ModelInfo.x0(:,1), ModelInfo.x0(:,2), 's','MarkerFaceColor','k', 'MarkerEdgeColor','k','MarkerSize',10);
    
%     buf{1} = sprintf('Low-fidelity training data (%d points)',size(ModelInfo.x1,1));
%     buf{2} = sprintf('High-fidelity training data (%d points)',size(ModelInfo.x2,1));
%     if time_dep == 1
%         buf{3} = sprintf('Anchor points (%d points)',size(ModelInfo.x0,1));
%     else
%         buf{3} = sprintf('Anchor points (%d points)',size(ModelInfo.x0,1));
%     end
%     
%     hl = legend(h,buf,'Location','southoutside'); legend boxoff
%     set(hl,'Interpreter','latex')
    
    %title('$\frac{|u_{2}(t,x) - \overline{u}_{2}(t,x)|}{||u_{2}||_{\infty}}$');
%     title('Error');
    colorbar
    if time_dep == 1
        xlabel('$t$')
        ylabel('$x$')
    else
        xlabel('$x_1$');
        ylabel('$x_2$');
    end
    axis square
    set(gca,'FontSize',48);
    set(gcf, 'Color', 'w');
    if save_plots  == 1
        export_fig error.png -r300
    end
    
    %     figure(5)
    %     clf
    %     hold
    % %     pcolor(Xplot,Yplot,corr(EH,EL));
    % %     shading interp
    % %     colorbar
    %     h(1) = plot(Exact_rhs_L, Exact_rhs_H,'bo','LineWidth',1);
    %     h(2) = plot(Exact_rhs_H, Exact_rhs_H,'k','LineWidth',1);
    %     xlabel('$f_{1}(x)$');
    %     ylabel('$f_{2}(x)$');
    %     title('$f_{1}(x),f_{2}(x)$ Cross-correlation');
    %     axis square
    %     set(gca,'FontSize',48);
    %     set(gca,'Xtick',[]);
    %     set(gca,'Ytick',[]);
    %     set(gcf, 'Color', 'w');
    %     if save_plots  == 1
    %         export_fig correlation.png -r300
    %     end
    
    %
    %     figure(3)
    %     clf
    %     hold
    %     contourf(Xplot,Yplot,abs(Exactplot - Predplot));
    %     scatter(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'r*');
    %     scatter(ModelInfo.x2(:,1), ModelInfo.x2(:,2), 'g*');
    %     colorbar
    %     title('Absolute error');
    %     axis square
    %
    %
    %     figure(4)
    %     clf
    %     hold
    %
    %     s1=surf(Xplot,Yplot,Exactplot_rhs);
    %     set(s1,'FaceColor',[0.0 0.2 0.0],'FaceAlpha',0.6);
    %     s2=surf(Xplot,Yplot,Predplot_rhs);
    %     set(s2,'FaceColor',[0 0 1],'FaceAlpha',0.6);
    %     legend('Exact RHS', 'Prediction RHS'); legend boxoff
    %     title('RHS prediction');
    %     axis square
    %
    %     figure(5)
    %     clf
    %     hold
    %     contourf(Xplot,Yplot,Varplot_rhs);
    %     scatter(ModelInfo.x1(:,1), ModelInfo.x1(:,2), 'r*');
    %     scatter(ModelInfo.x2(:,1), ModelInfo.x2(:,2), 'g*');
    %     title('RHS variance');
    %     colorbar
    %     axis square
    
end


