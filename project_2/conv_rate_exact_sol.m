%% convergence rate when known exact sol
clc;
close all;
clear;
%% compute error for ex. 1.1(b)
h = [0.04 0.02 0.01 0.005 0.0025];% 0.00125];
e = zeros(length(h), 1);
for i = 1:length(h)
    [cc, Xh, uh, e(i)] = main(h(i));
end
%% plot
close all;
x0=10;
y0=10;
width=500;
height=500;
hh = linspace(h(end),h(1));
figure;
loglog(h,e, '.', hh, hh.*7., 'k-','markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
legend('Error w.r.t. exact solution', 'f(dx)=dx');
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'FontSize',16);
xlim([1.e-3 5.e-2]);
ylim([1.e-2 5.e-1]);



