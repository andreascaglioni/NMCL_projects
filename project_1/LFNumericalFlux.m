function [FExt] = LFNumericalFlux(Uh, f, g)
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

    uh = Uh(:,2)./Uh(:,1);
    lMax =max((abs(uh)+(g*Uh(:,1)).^0.5));	% maximum wave speed 
    FExt = 0.5*(f(Uh(2:end,1),Uh(2:end,2))+f(Uh(1:end-1,1),Uh(1:end-1,2))...
          -0.5*lMax*(Uh(2:end,:)-Uh(1:end-1,:)));
end

