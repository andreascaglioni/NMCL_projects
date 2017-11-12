clc;
close all;
clear;
%% performance
dx = [0.01; 0.005; 0.0025; 0.00125; 0.000625];

 t1 = [0.021078 0.023511 0.020414 0.020510 0.021480;
       0.056424 0.050737 0.050780 0.051324 0.050160;
       0.185856 0.183104 0.197331 0.199605 0.195946;
       0.537020 0.538235 0.537102 0.523473 0.532027;
       1.709848 1.789554 1.766340 1.750415 1.797214];

 t2 = [0.048333 0.033018 0.036200 0.036147 0.039900;
     0.075271 0.075834 0.075253 0.077104 0.075289;
     0.277125 0.275365 0.279388 0.278109 0.281152;
     0.810744 0.809025 0.802739 0.807690 0.807626
     2.631252 2.615852 2.571127 2.547544 2.606592 ];
 
 t1Aver = zeros(length(dx),1);
 t2Aver = zeros(length(dx),1);
 for i = 1: size(dx)
     t1Aver(i) = sum(t1(i,:));
     t2Aver(i) = sum(t2(i,:));
 end
%% compute order of growth
compute_rate  = @(err, h)(log(err(1:end-1)./err(2:end)) ./ log(h(1:end-1)./h(2:end)));
p1 = compute_rate(t1Aver, dx);
averp1 = sum(p1)/length(p1);
p2 = compute_rate(t2Aver, dx);
averp2 = sum(p2)/length(p2);

%% plot
loglog(dx, t1Aver, '.', dx, t2Aver, '.', dx, dx.^(averp1)*1.e-4, 'k-',dx, dx.^(averp2)*1.e-4, 'k-', 'markersize', 15);
legend('first IC', 'second IC', num2str(averp1), num2str(averp2));
xlabel 'dx';
ylabel 'Average time';
grid on;