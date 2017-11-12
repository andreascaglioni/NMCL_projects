clc;
close all;
clear;
%% conpute error wrt referecen solution and plot
% refeence solution
dxRef = 0.000625;
[ccRef, hhRef, mhRef]=ex_1_2_a(dxRef);
%mesh spacings
dx = [0.01; 0.005; 0.0025; 0.00125];
err = zeros(length(dx),1);
for i = 1:length(dx)
    [cc, hh, mh]=ex_1_2_a(dx(i));
    % interpolate solution in reference mesh
    hh = repmat(hh, [dx(i)/dxRef 1]);
    hh = hh(:);
    mh = repmat(mh, [dx(i)/dxRef 1]);
    mh = mh(:);
    err(i) = dx(i)*(sum(abs(hh-hhRef))+sum(abs(mh-mhRef)));
end
%% plot
loglog(dx, err, '.',dx,dx.^(1)*5.e2, 'k-', 'markersize', 15);
grid on;
xlabel 'dx';
ylabel 'Error';
legend('Error wrt reference solution', 'f(dx)=dx');
