function [FExt] = LFNumericalFlux(u, v, f1, f2, g, maxVel)
% Computes Lax-Friedrichs numerical flux of finite volumes method applied to shallow water (1d). 
% data: 
%     hh    vecor of first unknown (height)
%     mh    vecor of second unknown (discharge)
%     f1    function handle of first component of the problem's flux
%     f2    function handle of secondcomponent of the problem's flux
%     g     real number (gravity acceleration)
% returns:
%   F1Ext   numerical flux for the first unknown (height) (computed also at interface with boundary conditions
%   F2Ext   numerical flux for the first unknown (discharge) (computed also at interface with boundary conditions
    
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

