%% convergence rate with reference numerical solution
clc;
close all;
clear;
%% conpute error wrt referecen solution and plot
% refence solution
dxRef = 0.000625;
[ccRef, hhRef, mhRef, uhRef, ~]=main(dxRef);
%mesh spacings
dx = [0.02; 0.01; 0.005; 0.0025; 0.00125];
err = zeros(length(dx),1);
for i = 1:length(dx)
    [cc, hh, mh, uh, ~]=main(dx(i));
    % interpolate solution in reference mesh
    hh = repmat(hh, [dx(i)/dxRef 1]);
    hh = hh(:);
    mh = repmat(mh, [dx(i)/dxRef 1]);
    mh = mh(:);
    err(i) = dx(i)*(sum(abs(hh-hhRef))+sum(abs(mh-mhRef)));
end
%% plot
x0=10;
y0=10;
width=500;
height=500;
close all;
figure;
loglog(dx, err, '.',dx,dx.^(1)*5.e2, 'k-', 'markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
legend('Error w.r.t. reference solution', 'f(dx)=dx');
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'FontSize',16);
xlim([1.e-3 3.e-2]);
ylim([3.e-1 15]);
    