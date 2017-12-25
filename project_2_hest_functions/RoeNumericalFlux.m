function [FExt] = RoeNumericalFlux(u, v)
% Computes Roe numerical flux of finite volumes method applied to shallow water (1d). 
% data: 
%     u,v       vectors argument of the flux
%     maxVel    maximum speed of the system
% returns:
%   FExt        vector containing the numerical flux for the unknowns by columns
%               (computed also at interface with boundary conditions)
    %% physical parameters
    g = 1;
    f1 = @(h,m) (m);
    f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
    %% compute
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

