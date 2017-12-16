function [FExt] = RoeNumericalFlux(u, v, f1, f2, g, maxVel)
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
    

    N = length(u)-1;
    % flux left
    hu = u(:,1);
    mu = u(:,2);
    Fu = [f1(hu,mu), f2(hu,mu)];
    % flux right
    hv = v(:,1);
    mv = v(:,2);
    Fv = [f1(hv,mv), f2(hv,mv)];
    
    %% change of variables
    z1 = 0.5*(hv.^0.5 + hu.^0.5);
    z1s = 0.5*(hv + hu);
    z2 = 0.5*(mv./((hv.^0.5)) + mu./((hu.^0.5)));
    w = z2./z1;
    %% Roe matrices
    ARoe = zeros(N+1,2,2); 
    for i = 1:N+1 
        A = [0., 1.; g*z1s(i)-w(i).^2., 2.*w(i)];
        [S,L] = eig(A);
        LRoe = abs(L);
        ARoe(i,:,:) = S*LRoe/S;
    end
    %% return the numerical flux
    FExt = 0.5*(Fu+Fv) ...
         - 0.5*([ARoe(:,1,1).*(hv-hu) + ARoe(:,1,2).*(mv-mu), ...
                 ARoe(:,2,1).*(hv-hu) + ARoe(:,2,2).*(mv-mu)]);
    
end

