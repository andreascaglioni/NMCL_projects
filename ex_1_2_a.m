function [cc, hh, mh, uh]=ex_1_2_a(dx)
% clc;
close all;
% clear;
%% choose BC
bcNumber = 0;   %DO NOT CHANGE
switch bcNumber
    case 0
        bcType = 'Periodic';
    case 1
        bcType = 'Open';
end
%% physical problem data
g = 1;

h0 = @(x) (1-0.1*sin(pi*x));
m0 = @(x) (0.*x);

% h0 = @(x) (1-0.2*sin(2*pi*x));
% m0 = @(x) (0.*x+0.5);

f1 = @(h,m) (m);
f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
%% domain and discretization
%dx = 0.01;  
xx = 0:dx:2;
cc = dx/2:dx:2-dx/2;
N = length(cc);
cc1=cc(1); ccN = cc(end);
switch bcNumber    %depending on type of BC extend cc (used to evaluate integral of S later)
    case 0
         ccExt = [ccN,cc,cc1];
    case 1
        ccExt = [cc1,cc,ccN];
end
CFL = 0.5;
T = 2.;
%% discrete IC (integrated)
hh0 = zeros(N, 1);      %stores discrete solution in a 1 dim vector (for every t)
mh0 = zeros(N, 1);
for j = 1:N
    hh0(j) = integral(h0,xx(j),xx(j+1),'AbsTol',1e-14)/dx;
    mh0(j) = integral(m0,xx(j),xx(j+1),'AbsTol',1e-14)/dx;
end
%% solution vector
hh = hh0;  %initialize solution to IC
mh = mh0;
uh = mh0./hh0;
tic;
%% solve
time = 0.;
k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
% i=1;
% figure;
while time < T
    if(time+k>T)
        time = T-k;
    end
    [hh, mh] = applyBC([hh;mh], bcType);     % enlarge the vectors by 2 each by applying BC
    uh = mh./hh;
    lMax = uh + (g*hh).^0.5;                 % local maximum wave speed (vector)
    
    % numerical flux
    F1Ext = 0.5*(f1(hh(2:end),mh(2:end)) + f1(hh(1:end-1),mh(1:end-1)))-0.5*lMax(1:end-1).*(hh(2:end) - hh(1:end-1));
    F2Ext = 0.5*(f2(hh(2:end),mh(2:end)) + f2(hh(1:end-1),mh(1:end-1)))-0.5*lMax(1:end-1).*(mh(2:end) - mh(1:end-1));
   % time advancement
    hh = hh(2:end-1) - k/dx*(F1Ext(2:end) - F1Ext(1:end-1));
    mh = mh(2:end-1) - k/dx*(F2Ext(2:end) - F2Ext(1:end-1));
    uh = mh./hh;
    
    time = time+k;
    %compute new timestep
    k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
    % store timestep in vector
    % kk(i)=k; i=i+1;
    % plot sl at every timestep
%     subplot(3,1,1)
%     plot(cc, hh);
%     grid on; title 'height';% ylim([0.9 1.1]);
%     subplot(3,1,2)
%     plot(cc, mh);
%     grid on; title 'discharge'; % ylim([-0.1 0.1]);
%     subplot(3,1,3)
%     plot(cc, uh); grid on; title 'speed';  %ylim([-0.1 0.1]);
%     pause(0.001);
end
toc
%% plot final solution
% figure
%     subplot(3,1,1)
%     plot(cc, hh);
%     grid on; title 'height'; %ylim([0.9 1.1]);
%     subplot(3,1,2)
%     plot(cc, mh);
%     grid on; title 'discharge'; % ylim([-0.03 0.03]);
%     subplot(3,1,3)
%     plot(cc, uh); grid on; title 'speed';  %ylim([-0.02 0.02]);
%% plot choices of k
% figure
% plot(kk, '.', 'markersize', 10);
% grid on;
% %% compute error at final time (in L^1)
% err = zeros(1,N);
% for j=1:N
%     hExaFinal =@(x) abs(hExa(x,T)-hh(j));
%     mExaFinal =@(x) abs(mExa(x,T)-mh(j));
%     err(j) = integral(hExaFinal,xx(j),xx(j+1),'absTol', 1e-14)...
%            + integral(mExaFinal,xx(j),xx(j+1),'absTol', 1e-14);
% end
% ERR = sum(err);
end
