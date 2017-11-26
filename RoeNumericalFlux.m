function [F1Ext, F2Ext] = RoeNumericalFlux(hh, mh, f1, f2, g)
% Computes Roe numerical flux of finite volumes method applied to shallow water (1d). 
% data: 
%     hh    vecor of first unknown (height)
%     mh    vecor of second unknown (discharge)
%     f1    function handle of first component of the problem's flux
%     f2    function handle of secondcomponent of the problem's flux
%     g     real number (gravity acceleration)
% returns:
%   F1Ext   numerical flux for the first unknown (height) (computed also at interface with boundary conditions
%   F2Ext   numerical flux for the first unknown (discharge) (computed also at interface with boundary conditions
    N = length(hh)-2;
    %% change of variables
    z1 = 0.5*(hh(2:end).^0.5 + hh(1:end-1).^0.5);
    z1s = 0.5*(hh(2:end) + hh(1:end-1));
    z2 = 0.5*(mh(2:end)./((hh(2:end).^0.5)) + mh(1:end-1)./((hh(1:end-1).^0.5)));
    w = z2./z1;
    %% Roe matrices
    ARoe = zeros(N+1,2,2); 
    for i = 1:N+1 
        A = [0., 1.; g*z1s(i)-w(i).^2., 2.*w(i)];
        [S,L] = eig(A);
        LRoe = abs(L);
        ARoe(i,:,:) = S*LRoe/S;
    end
    %% numerical flux
    F1Ext = 0.5*(f1(hh(2:end),mh(2:end)) + f1(hh(1:end-1),mh(1:end-1))) - 0.5*(ARoe(:,1,1).*(hh(2:end)-hh(1:end-1)) + ARoe(:,1,2).*(mh(2:end)-mh(1:end-1)));
    F2Ext = 0.5*(f2(hh(2:end),mh(2:end)) + f2(hh(1:end-1),mh(1:end-1))) - 0.5*(ARoe(:,2,1).*(hh(2:end)-hh(1:end-1)) + ARoe(:,2,2).*(mh(2:end)-mh(1:end-1)));
end

