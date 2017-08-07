clear all
clc
close all

set(0,'defaulttextinterpreter','latex')

eu = load('error_u.txt');
ef = load('error_f.txt');
hyps = load('hyps.txt');

Np = 80;
idx = [1 14 36 61];

ii = 4:1:Np+3;


color1 = [0 90 178]/255;
color2 = [208 61 3]/255;

figure(1)
clf
hold on
% plot(ii, eu_imse(1:101));
% plot(ii, ef_imse(1:101));
plot(ii(idx(1))*ones(10,1), linspace(1e-10,ef(idx(1)),10), '--k');
plot(ii(idx(1)), eu(idx(1)), 'o', 'MarkerSize', 14, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'LineWidth',2);
plot(ii(idx(1)), ef(idx(1)), 'o', 'MarkerSize', 14, 'MarkerEdgeColor', color1, 'LineWidth',2);

plot(ii(idx(2))*ones(10,1), linspace(1e-10,ef(idx(2)),10), '--k');
plot(ii(idx(2)), eu(idx(2)), 'o', 'MarkerSize', 14, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'LineWidth',2);
plot(ii(idx(2)), ef(idx(2)), 'o', 'MarkerSize', 14, 'MarkerEdgeColor', color1, 'LineWidth',2);

plot(ii(idx(3))*ones(10,1), linspace(1e-10,ef(idx(3)),10), '--k');
plot(ii(idx(3)), eu(idx(3)), 'o', 'MarkerSize', 14, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'LineWidth',2);
plot(ii(idx(3)), ef(idx(3)), 'o', 'MarkerSize', 14, 'MarkerEdgeColor', color1, 'LineWidth',2);

plot(ii(idx(4))*ones(10,1), linspace(1e-10,ef(idx(4)),10), '--k');
plot(ii(idx(4)), eu(idx(4)), 'o', 'MarkerSize', 14, 'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'LineWidth',2);
plot(ii(idx(4)), ef(idx(4)), 'o', 'MarkerSize', 14, 'MarkerEdgeColor', color1, 'LineWidth',2);

h(1) = plot(ii, eu(1:Np), 'Color', color2, 'LineWidth',3);
h(2) = plot(ii, ef(1:Np), '--', 'Color', color1, 'LineWidth',3);
xlabel('$n$')
ylabel('Relative error')
hl = legend(h, 'Solution', 'Forcing','Location','southoutside');
set(hl,'Interpreter','latex')



set(gca,'FontSize',24);
set(gcf, 'Color', 'w');
% set(gca,'YTick',1e-8:1e-2:1);

% legend('eu_imse','ef_imse','eu','ef');
set(gca,'YScale','log');
% 
% figure(2)
% clf
% hold all
% plot(ii, hyps(:,2));
% plot(ii, hyps(:,3));
% 
% figure(3)
% clf
% hold all
% plot(ii, hyps(:,1));
% set(gca,'YScale','log');
