function [ rhs ] = eval_rhs( Xh, f1, f2, g, S1Prim, S2Prim, xx, dx, time, maxVel, BCNumber, FluxNumber, SlopeNumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    N = size(Xh,1);
    % apply BC
    Xh = applyBC(Xh, BCNumber, 2);      % N+4
    % compute right and left finite differencs
    dp = Xh(3:N+4,:) - Xh(2:N+3,:);   % N+2
    dm = Xh(2:N+3,:) - Xh(1:N+2,:); 
    %compute slopes
    slope(:,1) = eval_slope(dm(:,1), dp(:,1), SlopeNumber); %N+2
    slope(:,2) = eval_slope(dm(:,2), dp(:,2), SlopeNumber);
    % compute current solutions at right and left cell interface
    Xhp = Xh(2:N+3,:) - 0.5*slope;      % N+2
    Xhm = Xh(2:N+3,:) + 0.5*slope;      % N+2
    % compute numerical flux
    Fh = eval_flux(Xhm(1:N+1,:), Xhp(2:N+2,:), FluxNumber, f1, f2, g,maxVel);     % N+1
    % compute forcing term
    uh = Xh(3:N+2,2)./Xh(3:N+2,1);  %N
    SAver = zeros(N,2);
    SAver(:,1) = (S1Prim(uh,xx(2:end)',time)-S1Prim(uh,xx(1:end-1)',time))/dx;
    SAver(:,2) = (S2Prim(uh,xx(2:end)',time)-S2Prim(uh,xx(1:end-1)',time))/dx;
    % return RHS
    rhs = - 1/dx*(Fh(2:N+1,:) - Fh(1:N,:)) + SAver;       % N
end

