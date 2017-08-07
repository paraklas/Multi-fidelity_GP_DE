clear all; Ntr = [0 60]; main;
clear all; Ntr = [180 0]; main;
clear all; Ntr = [180 60]; main;
clear all; Ntr = [360 120]; main;
clear all; Ntr = [180*3 180]; main;
close all
clear all
%clc

set(0,'defaulttextinterpreter','latex')
 
c1 = load('Kpred_d5_n1-0_n2-60.txt');
c2 = load('Kpred_d5_n1-180_n2-0.txt');
c3 = load('Kpred_d5_n1-180_n2-60.txt');
c4 = load('Kpred_d5_n1-360_n2-120.txt');
c5 = load('Exact_d5_n1-360_n2-120.txt');

% color5 =[27,158,119]/255;
% color4 =[217,95,2]/255;
% color3 =[117,112,179]/255;
% color2 =[231,41,138]/255;
% color1 =[102,166,30]/255;

% color4 =[27,158,119]/255;
% color3 =[217,95,2]/255;
% color2 =[117,112,179]/255;
% color1 =[231,41,138]/255;

% color4 = [106,81,163]/255;
% color3 = [158,154,200]/255;
% color2 = [203,201,226]/255;
% color1 = [242,240,247]/255;

% color1 =[200,200,200]/255;
% color2 =[150,150,150]/255;
% color3 =[80,80,80]/255;
% color4 =[40,40,40]/255;

% color1 =[255,255,178]/255;
% color2 =[254,204,92]/255;
% color3 =[253,141,60]/255;
% color4 =[227,26,28]/255;

color1 = [255,255,178]/255;
color2 = [254,204,92]/255;
color3 = [253,141,60]/255;
color4 = [240,59,32]/255;
color5 = [189,0,38]/255;

% color5 = [228,26,28]/255;
% color4 = [55,126,184]/255;
% color3 = [77,175,74]/255;
% color2 = [152,78,163]/255;
% color1 = [255,127,0]/255;



hold
s5 = histogram(c5,35,'Normalization','probability','FaceColor',color5,'FaceAlpha',1);
s4 = histogram(c4,35,'Normalization','probability','FaceColor',color4,'FaceAlpha',1);
s3 = histogram(c3,35,'Normalization','probability','FaceColor',color3,'FaceAlpha',1);
s2 = histogram(c2,35,'Normalization','probability','FaceColor',color2,'FaceAlpha',1);
s1 = histogram(c1,35,'Normalization','probability','FaceColor',color1,'FaceAlpha',1);




% 
% color4 = [106,81,163]/255;
% color3 = [158,154,200]/255;
% color2 = [203,201,226]/255;
% color1 = [242,240,247]/255;
% 
% color4 =[27,158,119]/255;
% color3 =[217,95,2]/255;
% color2 =[117,112,179]/255;
% color1 =[231,41,138]/255;




% color1 =[200,200,200]/255;
% color2 =[150,150,150]/255;
% color3 =[80,80,80]/255;
% color4 =[40,40,40]/255;
% 
% face_alpha = 0.5;
% 
% figure(1);
% clf
% hold all
% 
% plot(x1,e1,'color',color1,'LineWidth',2);
% s1 = stem(x1,e1,'color',color1);
% set(s1,'Marker','None');
% % set(s1, 'FaceColor',color1, 'FaceAlpha', 0.4);
% 
% plot(x2,e2,'color',color2,'LineWidth',2);
% s2 = stem(x2,e2,'color',color2);
% set(s2,'Marker','None');
% 
% % set(s2, 'FaceColor',color2, 'FaceAlpha', 0.4);
% 
% 
% plot(x3,e3,'color',color3,'LineWidth',2);
% s3 = stem(x3,e3,'color',color3);
% set(s3,'Marker','None');
% 
% % set(s3, 'FaceColor',color3, 'FaceAlpha', 0.4);
% 
% plot(x4,e4,'color',color4,'LineWidth',2);
% s4 = stem(x4,e4,'color',color4);
% set(s4,'Marker','None');

% set(s4, 'FaceColor',color4, 'FaceAlpha', 0.4);

% plot(x5,e5,'color',color5,'LineWidth',2);
% s5 = stem(x5,e5,'color',color5);
% set(s5,'Marker','None');

% set(s5, 'FaceColor',color5, 'FaceAlpha', 0.4);

% s1 = histogram(c1,25);
% s2 = histogram(c2,25);
% s3 = histogram(c3,25);
% s4 = histogram(c4,25);

buf = {'$n_1=0,n_2=60$', '$n_1=180,n_2=0$', '$n_1=180,n_2=60$','$n_1=360,n_2=120$','Exact'};
hl = legend([s1,s2,s3,s4,s5],buf,'Location','southoutside');
set(hl,'Interpreter','latex');


xlabel('$u_{2}(x), \overline{u}_{2}(x)$')
ylabel('frequency')
set(gca,'FontSize',24);
set(gcf, 'Color', 'w');