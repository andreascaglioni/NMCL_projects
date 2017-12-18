%% plots 3 numerical solutions
clc;
close all;
clear;
%% compute
[cc001, Xh001, uh001, ~] = main(0.01);
[cc0005, Xh0005, uh0005, ~] = main(0.005);
[ccRef, XhRef, uhRef, ~] = main(0.00125);
%% exact solution only for 1.1
u=0.25;
h0 = @(x) (1.+0.5*sin(pi*x));
hExa = @(x,t) h0(x-2);
mExa = @(x,t) u*hExa(x,2);
close all;
%% plots with OPTIMAL PLOT PARAMETERS
x0=10;
y0=10;
width=550;
height=700;
figure
subplot(2,1,1);
plot(ccRef(1:end-1), XhRef(:,1), cc001(1:end-1), Xh001(:,1), cc0005(1:end-1), Xh0005(:,1));
% hold on;
% plot(ccRef(1:end-1), hExa(ccRef(1:end-1)));
set(gcf,'units','points','position',[x0,y0,width,height]);
legend('dx=0.00125', 'dx=0.01', 'dx=0.005');%, 'Exact solution');
grid on;
title 'Height';
set(gca,'FontSize',16);
subplot(2,1,2);
plot(ccRef(1:end-1), XhRef(:,2), cc001(1:end-1), Xh001(:,2), cc0005(1:end-1), Xh0005(:,2));
% hold on;
% plot(ccRef(1:end-1), mExa(ccRef(1:end-1)));
legend( 'dx=0.00125', 'dx=0.01', 'dx=0.005');%, 'Exact solution');
grid on;
title 'Discharge';
set(gca,'FontSize',16);
% subplot(3,1,3);
% plot(ccRef(1:end-1), uhRef, cc001(1:end-1), uh001, cc0005(1:end-1), uh0005);
% legend('dx=0.00125', 'dx=0.01', 'dx=0.005');
% grid on;
% title 'Speed'