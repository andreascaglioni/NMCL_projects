function [F1Ext, F2Ext] = LFNumericalFlux(hh, mh, f1, f2, g)
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
    uh = mh./hh;
    lMax =max((abs(uh(1:end))+(g*hh(1:end)).^0.5));	% maximum wave speed (vector)
    F1Ext = 0.5*(f1(hh(2:end),mh(2:end)) + f1(hh(1:end-1),mh(1:end-1)))-0.5*lMax*(hh(2:end) - hh(1:end-1));
    F2Ext = 0.5*(f2(hh(2:end),mh(2:end)) + f2(hh(1:end-1),mh(1:end-1)))-0.5*lMax*(mh(2:end) - mh(1:end-1));
end

