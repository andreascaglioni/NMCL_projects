%% plots 3 numerical solutions
clc;
close all;
clear;
%% mesh size

%hh = [0.01 0.005 0.0025];
%hh = [0.04 0.02 0.01];
%hh = [0.08 0.04 0.02];
hh = [ 0.04 0.02 0.0025];
%% compute
[xx1, Xh1, uh1, ~] = main(hh(1));
[xx2, Xh2, uh2, ~] = main(hh(2));
[xx3, Xh3, uh3, ~] = main(hh(3));
%% exact solution only for 1.1
% u=0.25;
% h0 = @(x) (1.+0.5*sin(pi*x));
% hExa = @(x,t) h0(x-2);
% mExa = @(x,t) u*hExa(x,2);
% close all;
%% plots with OPTIMAL PLOT PARAMETERS
x0=10;
y0=10;
width=500;
height=700;
xxExa = 0:0.001:2;
figure
subplot(2,1,1);
plot( xx1, Xh1(:,1), xx2, Xh2(:,1), xx3, Xh3(:,1));
% hold on;
% plot(xxExa, hExa(xxExa));
set(gcf,'units','points','position',[x0,y0,width,height]);
legend( ['dx=' num2str(hh(1))], ['dx=' num2str(hh(2))],['dx=' num2str(hh(3))], 'Exact');
grid on;
title 'Height';
set(gca,'FontSize',16);
subplot(2,1,2);
plot( xx1, Xh1(:,2), xx2, Xh2(:,2), xx3, Xh3(:,2));
% hold on;
% plot(xxExa, mExa(xxExa));
legend( ['dx=' num2str(hh(1))], ['dx=' num2str(hh(2))],['dx=' num2str(hh(3))], 'Exact');
grid on;
title 'Discharge';
set(gca,'FontSize',16);
