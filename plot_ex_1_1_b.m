clc;
close all;
clear;
%% plot error in loglog
h = [0.02 0.01 0.005 0.0025 0.00125]
e = [   0.089257748248321 0.050078313615086 0.025252933916484 0.012513247992788    0.006072887959174]
hh = linspace(0.00125,0.02);
loglog(h,e, '.', hh, hh*2.5, 'k-','markersize', 15);
legend('data', 'h');
grid on;
xlabel 'h';
ylabel ' error';
xlim ([0., 0.025]);
ylim([0, 0.2]);
