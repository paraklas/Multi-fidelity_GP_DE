save_plots = 0;
close all
[X, Y] = meshgrid(linspace(min(nlow),max(nlow),80),linspace(min(nhigh),max(nhigh),80));
idx = (pp(:,2)==2);
Eplt = griddata(pp(~idx,1), pp(~idx,2), error_rhs(~idx), X, Y,'linear');
figure('units','normalized','outerposition',[0 0 1 1])
clf
hold
pcolor(X, Y, real(log10(Eplt)));
title('$\log_{10}(\frac{||u_2 - \overline{u}_2||_2}{||u_2||_2})$')
shading interp
colorbar
colormap cool

xlabel('$n_1$')
ylabel('$n_2$')
xlim([min(nlow) max(nlow)]);
ylim([min(nhigh+1) max(nhigh)]);
axis square
set(gca,'FontSize',24);
set(gcf, 'Color', 'w');
set(gca,'XTick',0:5:30)
set(gca,'YTick',4:1:7)

if save_plots  == 1
    export_fig convergence_rhs.png -r300
end

[X, Y] = meshgrid(linspace(min(nlow),max(nlow),80),linspace(min(nhigh),max(nhigh),80));
idx = (pp(:,2)==2);
Eplt = griddata(pp(~idx,1), pp(~idx,2), error_u(~idx), X, Y,'linear');
figure('units','normalized','outerposition',[0 0 1 1])
clf
hold
pcolor(X, Y, real(log10(Eplt)));
title('$\log_{10}(\frac{||f_2 - \overline{f}_2||_2}{||f_2||_2})$')
shading interp
colorbar
colormap cool

xlabel('$n_1$')
ylabel('$n_2$')
xlim([min(nlow) max(nlow)]);
ylim([min(nhigh+1) max(nhigh)]);
axis square
set(gca,'FontSize',24);
set(gcf, 'Color', 'w');
set(gca,'XTick',0:5:30)
set(gca,'YTick',4:1:7)

if save_plots  == 1
    export_fig convergence_u.png -r300
end