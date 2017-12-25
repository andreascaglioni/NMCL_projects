function [numflux] = EulerLF(u,v,maxvel)
% function [numflux] = EulerLF(u,v,gamma,lambda,maxvel);
% Purpose: Evaluate global Lax Friedrich numerical flux for 
% the Euler equations

%% physical parameters
    g = 1;
    f1 = @(h,m) (m);
    f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
    %% compute
    % flux left
    hu = u(:,1);
    mu = u(:,2);
    fu = [f1(hu,mu), f2(hu,mu)];
    % flux right
    hv = v(:,1);
    mv = v(:,2);
    fv = [f1(hv,mv), f2(hv,mv)];
% Evaluate numerical flux
numflux = (fu+fv)/2 - maxvel/2*(v-u);
return