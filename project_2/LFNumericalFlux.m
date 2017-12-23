function [FExt] = LFNumericalFlux(u, v, maxVel)
% Computes Lax-Friedrichs numerical flux of finite volumes method applied to shallow water (1d). 
% data: 
%     u,v       vectors argument of the flux
%     maxVel    maximum speed of the system
% returns:
%     FExt      vector containing the numerical flux for the unknowns by columns
%               (computed also at interface with boundary conditions)
    %% physical parameters
    g = 1;
    f1 = @(h,m) (m);
    f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
    %% compute
    % flux left
    hu = u(:,1);
    mu = u(:,2);
    Fu = [f1(hu,mu), f2(hu,mu)];
    % flux right
    hv = v(:,1);
    mv = v(:,2);
    Fv = [f1(hv,mv), f2(hv,mv)];
    % return flux
    FExt = 0.5*(Fu+Fv) - 0.5*maxVel*(v-u);
end

