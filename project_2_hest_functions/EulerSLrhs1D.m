function [dq] = EulerSLrhs1D(q,x,h,time,maxvel, BCNumber, FluxNumber, SlopeNumber, SourceNumber, c)
% function [dq] = EulerSLrhs1D(x,q,gamma,h,k,maxvel);
% Purpose: Evaluate right hand side for the Euler equations using 
% slope-limited method with limiting on characteristic variables
N = length(x);
dq = zeros(N,2);
qL = zeros(N+2,2); qR = zeros(N+2,2);
dup = zeros(N+2,2); dum = zeros(N+2,2); 
duL = zeros(N+2,2);

switch BCNumber
    case 0
        BCString = 'P';
    case 1
        BCString = 'N';
end
% Extend data and assign boundary conditions
[xe,re] = extend(x,q(:,1),h,2,BCString,0,BCString,0);
[xe,me] = extend(x,q(:,2),h,2,BCString,0,BCString,0);

% Extract variables
qe = [re me];

% Compute left and right differences, evaluate slopes and interface values
dup = qe(3:N+4,:) - qe(2:N+3,:);
dum = qe(2:N+3,:) - qe(1:N+2,:);
duL(:,1) = SlopeLimit(dup(:,1), dum(:,1),SlopeNumber);
duL(:,2) = SlopeLimit(dup(:,2), dum(:,2),SlopeNumber);
qL = qe(2:N+3,:) - duL/2;
qR = qe(2:N+3,:) + duL/2;
 
%% source term
g = 1;
switch SourceNumber 
    case 0
        %S1 = @(v, x,t) (pi/2*(v-1).*cos(pi*(x-t)));
        %S2 = @(v, x,t) (pi/2*cos(pi*(x-t)).*(-v+v.^2+g*h0(x-t)));
        S1Prim = @(v, x,t) (0.5*(v-1).*sin(pi*(x-t)));
        S2Prim = @(v, x,t) ( 0.5*(g-v+v.^2).*sin(pi*(x-t)) ...
                           + 0.25*g*pi* (0.5/pi*sin(pi*(x-t)).^2));
    case 1
        %S1 = @(v, x,t) 0*v+0*x+0*t;
        %S2 = @(v, x,t)  0*v+0*x+0*t;
        S1Prim = @(v, x,t)  0*v+0*x+0*t;
        S2Prim = @(v, x,t)  0*v+0*x+0*t;
end
% compute forcing term
    uh = q(:,2)./q(:,1);
    SAver = [(S1Prim(uh,c(2:end)',time)-S1Prim(uh,c(1:end-1)',time))/h, ...
             (S2Prim(uh,c(2:end)',time)-S2Prim(uh,c(1:end-1)',time))/h]; %N
%%  
         
% Evaluate right hand side using numerical flux
switch FluxNumber
    case 0
        dq = - (EulerLF(qR(2:N+1,:),qL(3:N+2,:),maxvel) - ...
                EulerLF(qR(1:N,:),qL(2:N+1,:),maxvel))/h + SAver;
    case 1
        dq = - (RoeNumericalFlux(qR(2:N+1,:),qL(3:N+2,:)) - ...
                RoeNumericalFlux(qR(1:N,:),qL(2:N+1,:)))/h + SAver;
end
return
