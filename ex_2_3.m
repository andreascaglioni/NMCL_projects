clc;
close all;
clear;
%% choose BC
bcNumber = 0;
switch bcNumber
    case 0
        bcType = 'Periodic';
    case 1
        bcType = 'Open';
end
%% physical problem data
g = 1;
u = 0.25;   %INITIAL velocity
h0 = @(x) (1+0.5*sin(pi*x));
m0 = @(x) (u*h0(x));
f1 = @(h,m) (m);
f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
S1 = @(v, x,t) (pi/2*(v-1).*cos(pi*(x-t)));
S2 = @(v, x,t) (pi/2*cos(pi*(x-t)).*(-v+v.^2+g*h0(x-t)));
S1Prim = @(v, x,t) (0.5*(v-1).*sin(pi*(x-t)));
S2Prim = @(v, x,t) ( 0.5*(g-v+v.^2).*sin(pi*(x-t)) ...
                   + 0.25*g*pi* (0.5/pi*sin(pi*(x-t)).^2));
%% exact solution
hExa = @(x,t) h0(x-t);
mExa = @(x,t) u*hExa(x,t);
%% domain and discretization
dx = 0.01;  
xx = 0:dx:2;
cc = dx/2:dx:2-dx/2;
N = length(cc);
cc1=cc(1); ccN = cc(end);
switch bcNumber        %depending on type of BC extend cc (used to evaluate integral of S later)
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
%% solution matrix
hh = hh0;  %initialize solution to IC
mh = mh0;
uh = mh0./hh0;
%% solve
k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
time = 0.;
i=1;
figure(1);
while time < T
    if(time+k>T)
        time = T-k;
    end
    [hh, mh] = applyBC([hh;mh], bcType);     % enlarge the vectors by 2 each by applying BC
    uh = mh./hh;
    %% computing Roe matrix and numerical flux
    z1 = 0.5*(hh(2:end).^0.5 + hh(1:end-1).^0.5);
    z1s = 0.5*(hh(2:end) + hh(1:end-1));
    z2 = 0.5*(mh(2:end)./((hh(2:end).^0.5)) + mh(1:end-1)./((hh(1:end-1).^0.5)));
    w = z2./z1;
    F1Ext = zeros(N+1,1);
    F2Ext = zeros(N+1,1);
    for i = 1:N+1
        lambdaMax = w(i)+(g*z1s(i))^0.5;
        lambdaMin = w(i)-(g*z1s(i))^0.5;
        vMax = 1/(1+(w(i)+(g*z1s(i))^0.5)^2)^0.5*[1; w(i)+(g*z1s(i))^0.5];
        vMin = 1/(1+(-w(i)+(g*z1s(i))^0.5)^2)^0.5*[1;w(i)-(g*z1s(i))^0.5];
        SRoe = [vMin, vMax];
        LRoe = [abs(lambdaMin), 0; 0, abs(lambdaMax)];
        ARoe = SRoe\LRoe*SRoe;        
        % numerical flux
        F1Ext(i) = 0.5*(f1(hh(i+1),mh(i+1)) + f1(hh(i),mh(i))) - 0.5*(ARoe(1,:)*[hh(i+1)-hh(i); mh(i+1)-mh(i)]);
        F2Ext(i) = 0.5*(f2(hh(i+1),mh(i+1)) + f2(hh(i),mh(i))) - 0.5*(ARoe(2,:)*[hh(i+1)-hh(i); mh(i+1)-mh(i)]);
    end
    % local average of the RHS
    S1Int = (S1Prim(uh(2:end), ccExt(2:end)', time) - S1Prim(uh(1:end-1), ccExt(1:end-1)', time))/dx;
    S2Int = (S2Prim(uh(2:end), ccExt(2:end)', time) - S2Prim(uh(1:end-1), ccExt(1:end-1)', time))/dx;
    % time advancement
    hh = hh(2:end-1) - k/dx*(F1Ext(2:end) - F1Ext(1:end-1)) + k*S1Int(1:end-1);
    mh = mh(2:end-1) - k/dx*(F2Ext(2:end) - F2Ext(1:end-1)) + k*S2Int(1:end-1);
    uh = mh./hh;
    time = time+k;
    %compute new timestep
    k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
    % store timestep in vector
    % kk(i)=k; i=i+1;
    % plot sl at every timestep
    subplot(3,1,1)
    plot(cc, hh, cc ,hExa(cc,time));
    grid on; title 'height';
    subplot(3,1,2)
    plot(cc, mh, cc ,mExa(cc,time));
    grid on; title 'discharge';
    subplot(3,1,3)
    plot(cc, uh);
    grid on; title 'speed';
%     pause(0.001);
%     press=waitforbuttonpress;
end
%% plot final solution
figure
    subplot(3,1,1)
    plot(cc, hh, cc ,hExa(cc,time));
    grid on; title 'height';
    subplot(3,1,2)
    plot(cc, mh, cc ,mExa(cc,time));
    grid on; title 'discharge';
    subplot(3,1,3)
    plot(cc, uh); grid on; title 'speed';
%% plot choices of k
% figure
% plot(kk, '.', 'markersize', 10);
% grid on;
%% compute error at final time (in L^1)
% err = zeros(1,N);
% for j=1:N
%     hExaFinal =@(x) abs(hExa(x,T)-hh(j));
%     mExaFinal =@(x) abs(mExa(x,T)-mh(j));
%     err(j) = integral(hExaFinal,xx(j),xx(j+1),'absTol', 1e-14)...
%            + integral(mExaFinal,xx(j),xx(j+1),'absTol', 1e-14);
% end
% ERR = sum(err);
