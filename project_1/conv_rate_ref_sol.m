%% convergence rate with reference numerical solution
clc;
close all;
clear;
%% conpute error wrt referecen solution and plot
% reference solution
dxRef =  0.000625;
[ccRef,UhRef, uhRef, ~]=main(dxRef);
%mesh spacings
dx = [0.04; 0.02; 0.01; 0.005; 0.0025];% 0.00125];
err = zeros(length(dx),1);
for i = 1:length(dx)
    [cc, Uh, uh, ~]=main(dx(i));
    % interpolate solution in reference mesh
    hhInt = repmat(Uh(:,1), [1, dx(i)/dxRef])';
    hhInt = hhInt(:);
    mhInt = repmat(Uh(:,2), [1, dx(i)/dxRef])';
    mhInt = mhInt(:);
    err(i) = dxRef*(sum(abs(hhInt-UhRef(:,1)))+sum(abs(mhInt-UhRef(:,2))));
end
%% plot
x0=10; y0=10; width=500; height=500;
close all;
figure;
loglog(dx, err, '.',dx,1.5*dx.^1, 'k-', 'markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
legend('Error w.r.t. reference solution', 'f(dx)=dx');
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'FontSize',16);
xlim([1.e-3 5.e-2]);
%%ylim([1.e-2 3.e-1]);
    