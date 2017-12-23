function [ rhs ] = eval_rhs( Xh, cc, dx, time, maxVel, BCNumber, FluxNumber, SlopeNumber, SourceNumber)
% Evaluates the right hand side of the equation
% Data:
%     Xh            vector containing the solution at current time/ RK step
%     cc            vector of cells interfaces
%     dx            mesh spacing
%     time          current time
%     maxVel        maximum velocity of the system
%     BCNumber      number to choose boundary condition
%     FluxNumber    number to choose the numerical flux
%     SlopeNumber   number to choose which slope to use
%     SourceNumber  number to choose the source term of the problem
% Returns:
%     rhs           vector containing the right hand side  
%% source term
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
%% computing
    % compute forcing term
    N = size(Xh,1);   % number of cells
    uh = Xh(:,2)./Xh(:,1);
    SAver = [(S1Prim(uh,cc(2:end)',time)-S1Prim(uh,cc(1:end-1)',time))/dx, ...
             (S2Prim(uh,cc(2:end)',time)-S2Prim(uh,cc(1:end-1)',time))/dx]; %N
    % apply BC
    Xh = applyBC(Xh, BCNumber, 2);      % N+4
    % compute right and left finite differencs
    dp = Xh(3:N+4,:) - Xh(2:N+3,:);   % N+2
    dm = Xh(2:N+3,:) - Xh(1:N+2,:); 
    %compute slopes
    slope = eval_slope(dp, dm, SlopeNumber); %N+2
    % compute current solutions at right and left cell interface
    Xhp = Xh(2:N+3,:) - 0.5*slope;      % N+2
    Xhm = Xh(2:N+3,:) + 0.5*slope;      % N+2
    % return RHS
    rhs = -1/dx*(eval_flux(Xhm(2:N+1,:), Xhp(3:N+2,:), FluxNumber, maxVel) ...
               - eval_flux(Xhm(1:N,:), Xhp(2:N+1,:), FluxNumber, maxVel)) + SAver;
end