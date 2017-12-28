%% convergence rate with reference numerical solution
clc;
close all;
clear;
%% conpute error wrt referecen solution and plot
% reference solution
dxRef =  0.00125;
[xxRef, XhRef,~, ~]=main(dxRef);
%XhRef = XhRef(1:end-1,:); 
hhRef = XhRef(:,1); mhRef=XhRef(:,2);
%mesh spacings
dx = [0.04; 0.02; 0.01; 0.005; 0.0025];% 0.00125];
err = zeros(length(dx),1);
for i = 1:length(dx)
    [xx, Xh, ~, ~]=main(dx(i));
    %Xh = Xh(1:end-1,:);
    % interpolate solution in reference mesh
    hhInt = interp1(xx, Xh(:,1), xxRef);
    mhInt = interp1(xx, Xh(:,2), xxRef);
    err(i) = dxRef*(sum(abs(hhRef-hhInt')) + sum(abs(mhRef-mhInt')));
end
%% plot
x0=10;
y0=10;
width=500;
height=500;
close all;
figure;
loglog(dx, err, '.',dx,6*dx.^1, 'k-', 'markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
legend('Error w.r.t. reference solution', 'f(dx)=dx');
set(gcf,'units','points','position',[x0,y0,width,height]);
set(gca,'FontSize',16);
% xlim([1.e-3 5.e-2]);
% ylim([1.e-2 3.e-1]);