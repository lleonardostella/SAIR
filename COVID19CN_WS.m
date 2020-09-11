clear; clc; clf;
delete(findall(0,'Type','figure'))

N = 1000;
m = 4;
beta = .35;

figure
h4 = WattsStrogatz(N,m,beta);
plot(h4,'NodeColor','k','EdgeAlpha',0.5);
title('Watts-Strogatz Complex Network','FontSize',16);
set(gca,'FontSize',14);

figure
colormap default
deg = degree(h4);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(h4,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.5)
title('Watts-Strogatz Graph - Node Degree by Color','FontSize',16);
colorbar

figure
colormap spring
deg = degree(h4);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(h4,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.5)
title('Watts-Strogatz Graph - Node Degree by Color','FontSize',16);
colorbar

figure
colormap bone
deg = degree(h4);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(h4,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.5)
title('Watts-Strogatz Graph - Node Degree by Color','FontSize',16);
colorbar

figure
colormap hot
deg = degree(h4);
nSizes = 2*sqrt(deg-min(deg)+0.2);
nColors = deg;
plot(h4,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.5)
title('Watts-Strogatz Graph - Node Degree by Colour','FontSize',16);
colorbar

figure
%histogram(degree(h2),'BinMethod','integers','FaceAlpha',0.9);
%hold on
%histogram(degree(h3),'BinMethod','integers','FaceAlpha',0.9);
histogram(degree(h4),'BinMethod','integers','facecolor',[0.2 0.5 0.9],'FaceAlpha',0.8);
title('Node Degree Distributions for Watts-Strogatz Model','FontSize',16);
xlabel('Node degree','FontSize',14);
ylabel('Number of nodes','FontSize',14);
set(gca,'FontSize',14);
legend('p = 1','Location','NorthWest');