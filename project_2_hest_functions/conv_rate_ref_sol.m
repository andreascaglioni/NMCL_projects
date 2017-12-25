%% convergence rate with reference numerical solution
% clc;
% close all;
% clear;
%% conpute error wrt referecen solution and plot
% reference solution
dxRef =  0.0025 ;%0.000625;
[ccRef, XhRef, uhRef, ~]=main(dxRef);
hhRef = XhRef(:,1); mhRef=XhRef(:,2);
%mesh spacings
dx = [0.04 0.02 0.01 0.005];
err = zeros(length(dx),1);
for i = 1:length(dx)
    [cc, Xh, uh, ~]=main(dx(i));
    hh = Xh(:,1); mh = Xh(:,2);
    % interpolate solution in reference mesh
    hhInt = repmat(hh, [1, dx(i)/dxRef])'; 
    hhInt = hhInt(:);
    mhInt = repmat(mh, [1, dx(i)/dxRef])';
    mhInt = mhInt(:);
    err(i) = dxRef*(sum(abs(hhInt-hhRef))+sum(abs(mhInt-mhRef)));
end
%% plot
x0=10;
y0=10;
width=500;
height=500;
% close all;
% figure;
% loglog(dx, 4*dx, 'k-');
% hold on;
loglog(dx, err, '.', 'markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
% legend('Error w.r.t. reference solution', 'f(dx)=dx');
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'FontSize',16);
% xlim([1.e-3 5.e-2]);
% ylim([1.e-2 3.e-1]);
hold on;
    