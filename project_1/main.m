function [cc, Uh, uh, ERR] = main(dx)
% the "main" function fo the assignment which executes the task and calls
% the other functions. The only data is the meshs pacing, howoever all
% other parametnes can be chose in the section "choose" (boundary 
% conditions, initial conditions, source term, exact solution, numerical 
% flux) (see comments) and " domain and discretization" (CFL number, final
% time).
% data:
%     dx      spacing for space discretization
% returns:    
%     cc      vector containig interface points of cells
%     hh      numerical solution for height at final time
%     mh      numerical solution for discharge at final time
%     uh      numerical solution for speed (ratio discharge/height) at 
%             final time
%     ERR     L1 error of the numerical solution wrt the analytical 
%        solution (n.b. relevant only if the analytical solution is available
%% choose
BCNumber = 0;       % 0:periodic;   1:open
ICNumber = 1;       % 0: ex. 1.1;   1,2: ex. 1.2(a) and 1.2(b);     3: ex. 1.4
SourceNumber = 1;   % 0: ex: 1.1;   1: no source
ExaNumber = 1;      % 0:ex. 1.1;    1: no exact sol. available (put it to 0)
FluxNumber = 0;     % 0:LF;         1:Roe
%% physical parameters
g = 1;
f = @(h,m) [m, m.^2./h + 0.5*g*h.^2];
%% IC
switch ICNumber
    case 0
        u = 0.25;   %INITIAL velocity
        h0 = @(x) (1.+0.5*sin(pi*x));
        m0 = @(x) (u*h0(x));
    case 1
        h0 = @(x) (1.-0.1*sin(pi*x));
        m0 = @(x) (0.*x);
    case 2
        h0 = @(x) (1.-0.2*sin(2*pi*x));
        m0 = @(x) (0*x+0.5);
    case 3
        h0 = @(x) (1.+0*x);
        m0 = @(x) (-1.5*(x<1.) -0.75*(x==1.));
end
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
%% exact solution
switch ExaNumber
    case 0
        hExa = @(x,t) h0(x-t);
        mExa = @(x,t) u*hExa(x,t);
    case 1
        hExa = @(x,t) 0*x+0*t;
        mExa = @(x,t) 0*x+0*t;
end
%% domain and discretization
cc = 0:dx:2;            %cells boundaries
xx = dx/2:dx:2-dx/2;    %cells' centerpoints
N = length(xx);         % number of cells
CFL = 0.5;
T = 2.;
%% discrete IC (integrated)
Uh0 = zeros(N, 2);      %stores averaged IC in a 1 dim vector (for every t)
for j = 1:N
    Uh0(j,1) = integral(h0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
    Uh0(j,2) = integral(m0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
end
%% solution vectors
Uh = Uh0;        %initialize solution to IC
uh = Uh(:,2)./Uh(:,1);
% plotSOlAtTIme(xx, hh, @(x)hExa(x,0.), mh, @(x)mExa(x,0.));
% press  = waitforbuttonpress;
%% solve
time = 0.;
% figure;
while time < T
    k = CFL*dx/max(abs(uh)+(g*Uh(:,1)).^0.5);
    if(time+k>T)
        time = T-k;
    end
    Uh = applyBC(Uh, BCNumber);   % enlarge the vectors by 2 each by applying BC
    uh = Uh(:,2)./Uh(:,1);        % also extend u
    % numerical flux
    switch FluxNumber
        case 0
            FExt = LFNumericalFlux(Uh, f, g);
        case 1
            FExt = RoeNumericalFlux(Uh, f, g);
    end
    % local average of the RHS
    S1Int = (S1Prim(uh(2:end-1),cc(2:end)',time)-S1Prim(uh(2:end-1),cc(1:end-1)',time))/dx;
    S2Int = (S2Prim(uh(2:end-1),cc(2:end)',time)-S2Prim(uh(2:end-1),cc(1:end-1)',time))/dx;
    % time advancement
    Uh = Uh(2:end-1,:) - k/dx*(FExt(2:end,:)-FExt(1:end-1,:)) + k*[S1Int, S2Int];
    uh = Uh(:,2)./Uh(:,1);
    time = time+k;
    
    %plot sol at every timestep
%   plotSolAtTIme(xx, Uh, @(x)hExa(x,time), @(x)mExa(x,time));
%   press = waitforbuttonpress;
%   pause(0.001);
end
%% plot final solution
figure
plotSolAtTIme(xx, Uh, @(x)hExa(x,T), @(x)mExa(x,T));
%% compute error at final time (in L^1)
ERR = computeL1Error(hExa, mExa, Uh, cc, T);
end
