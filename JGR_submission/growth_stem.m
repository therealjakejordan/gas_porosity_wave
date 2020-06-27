clear all
% close all,
clc

addpath('./Colonel/Hamburger_Helper')
addpath('./Colonel/Smooth_Operators')

set(groot,'defaultAxesTickLabelinterpreter','latex')
load('./test_processed.mat')
freq = ff;
iz   = indexZ;
fnt = 28;
%% colors
yale1 = [0 53 107]/255;
yale2 = [189 83 25]/255;
yale3 = [99 170 255]/255;
yale4 = [95 113 45]/255;
yale5 = [74 74 74]/255;
yale6 = [40 109 192]/255;

%%
figure('units','normalized','outerposition',[0.0 0.0 0.6 0.6])
color = orangeblue(128);
 
stem3(FM(1:40:end,freq(2:end)),XCf(1:40:end,freq(2:end)),PmcM(iz(1:40:end),freq(2:end)),...
    'Color',[150 151 152]/255), hold on
stem3(FM(1:40:end,freq(3)),XCf(1:40:end,freq(3)),PmcM(iz(1:40:end),freq(3)),...
    'Color',yale3,'MarkerFaceColor',yale3,'MarkerEdgeColor',yale3) %3 2 5
stem3(FM(1:40:end,freq(5)),XCf(1:40:end,freq(5)),PmcM(iz(1:40:end),freq(5)),...
    'Color',yale2,'MarkerFaceColor',yale2,'MarkerEdgeColor',yale2)
stem3(FM(1:40:end,freq(6)),XCf(1:40:end,freq(6)),PmcM(iz(1:40:end),freq(6)),...
    'Color',yale5,'MarkerFaceColor',yale5,'MarkerEdgeColor',yale5)
stem3(FM(1:40:end,freq(8)),XCf(1:40:end,freq(8)),PmcM(iz(1:40:end),freq(8)),...
    'Color',yale4,'MarkerFaceColor',yale4,'MarkerEdgeColor',yale4)
set(gca,'TickLabelInterpreter','latex','fontsize',18)

xlim([0 102]/char)
ylim([0 113])
view([45 45])
xlabel('$f$','interpreter','latex','fontsize',24)

ylabel('$z$','interpreter','latex','fontsize',24)
zlabel('$\Pi$','interpreter','latex','fontsize',24)
zlim([0 2]*1e-5)
 
fig1=figure(1);
fig1.Renderer='Painters';