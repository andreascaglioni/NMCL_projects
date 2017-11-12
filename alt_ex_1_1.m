clc;
close all;
clear;

%% choose bc
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
N = length(xx);
CFL = 0.5;
T = 2.;
%% discrete IC (integrated)
hh0 = h0(xx)';      %stores discrete solution in a 1 dim vector (for every t)
mh0 = m0(xx)';

%% solution matrixh
hh = hh0;  %initialize solution to IC
mh = mh0;
uh = u*ones(N,1);
k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
%% solve
time = 0.;
i=1;
figure(1);
while time < T
    time
    if(time+k>T)
        time = T-k;
    end
    
    [hh, mh] = applyBC([hh;mh], bcType);          % enlarge the vectors by 2 each by applying BC

    MaxV = u - (g*hh).^0.5;                 % local maximum wave speed (vector)
    
    F1Ext = 0.5*(f1(hh(2:end),mh(2:end)) - f1(hh(1:end-1),mh(1:end-1))) - 0.5*MaxV(2:end).*(hh(2:end) - hh(1:end-1));
    F2Ext = 0.5*(f2(hh(2:end),mh(2:end)) - f2(hh(1:end-1),mh(1:end-1))) - 0.5*MaxV(2:end).*(mh(2:end) - mh(1:end-1));
    
    hh = hh(2:end-1) - k/dx*(F1Ext(2:end) - F1Ext(1:end-1)) + k*S1(uh,xx',time);
    mh = mh(2:end-1) - k/dx*(F2Ext(2:end) - F2Ext(1:end-1)) + k*S2(uh,xx',time);
    
    uh = mh./hh;
    
    time = time+k;
    
    k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
%   kk(i)=k; i=i+1;
    subplot(3,1,1)
    plot(xx, hh, xx ,hExa(xx,time));
    grid on;
    title 'height'
    subplot(3,1,2)
    plot(xx, mh, xx ,mExa(xx,time));
    grid on;
    title 'discharge'
    subplot(3,1,3)
    plot(xx, uh);
    grid on;
    title 'speed'
    pause(0.001);
    %if(time > 0.6)
%         press = waitforbuttonpress;
    %end
end
% figure
%     subplot(2,1,1)
%     plot(cc, hh);
%     grid on;
%     ylim([-2 2]);
%     title 'height'
%     subplot(2,1,2)
%     plot(cc, mh);
%     ylim([-2 2]);
%     grid on;
%     title 'discharge'
% figure
% plot(kk, '.', 'markersize', 10);
