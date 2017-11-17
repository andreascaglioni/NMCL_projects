function [cc, hh, mh, uh, ERR]=main(dx)
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
BCNumber = 0;   % 0:periodic;   1:open
ICNumber = 0;   % 0: ex. 1.1;   1,2: ex. 1.2(a) and 1.2(b);     3: ex. 1.4
SourceNumber = 0;  % 0: ex: 1.1;   1: no source
ExaNumber = 0;  % 0:ex. 1.1;   1: no exact sol available (put it to 0)
FluxNumber = 1; % 0:LF;     1:Roe
%% BC string
switch BCNumber
    case 0
        bcType = 'Periodic';
    case 1
        bcType = 'Open';
end
%% physical parameters
g = 1;
f1 = @(h,m) (m);
f2 = @(h,m) (m.^2./h + 0.5*g*h.^2);
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
        m0 = @(x) (-1.5*(x<1));
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
cc = 0:dx:2;    %cells boundaries
xx = dx/2:dx:2-dx/2;    %cells' centerpoints
N = length(xx); % number of cells
CFL = 0.5;
T = 2.;
%% discrete IC (integrated)
hh0 = zeros(N, 1);      %stores discrete solution in a 1 dim vector (for every t)
mh0 = zeros(N, 1);
for j = 1:N
    hh0(j) = integral(h0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
    mh0(j) = integral(m0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
end
%% solution vectors
hh = hh0;  %initialize solution to IC
mh = mh0;
uh = mh0./hh0;
% plotSOlAtTIme(xx, hh, @(x)hExa(x,0.), mh, @(x)mExa(x,0.));
% press  = waitforbuttonpress;
%% solve
k = CFL*dx/max(abs(uh)+(g*hh).^0.5);    %determine timestep
time = 0.;
% figure;
while time < T
    if(time+k>T)
        time = T-k;
    end
    [hh, mh] = applyBC(hh, mh, bcType);     % enlarge the vectors by 2 each by applying BC
    uh = mh./hh;                            %also extend u
    % numerical flux
    switch FluxNumber
        case 0
            [F1Ext, F2Ext] = LFNumericalFlux(hh, mh, f1, f2, g);
        case 1
            [F1Ext, F2Ext] = RoeNumericalFlux(hh, mh, f1, f2, g);
    end
    % local average of the RHS
    S1Int = (S1Prim(uh(2:end-1), cc(2:end)', time) - S1Prim(uh(2:end-1), cc(1:end-1)', time))/dx;
    S2Int = (S2Prim(uh(2:end-1), cc(2:end)', time) - S2Prim(uh(2:end-1), cc(1:end-1)', time))/dx;
    % time advancement
    hh = hh(2:end-1) - k/dx*(F1Ext(2:end) - F1Ext(1:end-1)) + k*S1Int;
    mh = mh(2:end-1) - k/dx*(F2Ext(2:end) - F2Ext(1:end-1)) + k*S2Int;
    uh = mh./hh;
    time = time+k;
    %compute new timestep
    k = CFL*dx/max(abs(uh)+(g*hh).^0.5);
    %plot sol at every timestep
%     plotSOlAtTIme(xx, hh, @(x)hExa(x,time), mh, @(x)mExa(x,time));
%     press = waitforbuttonpress;
%     pause(0.001);
end
%% plot final solution
figure
plotSOlAtTIme(xx, hh, @(x)hExa(x,T), mh, @(x)mExa(x,T));
%% compute error at final time (in L^1)
ERR = computeL1Error(hExa, mExa, hh, mh, cc, T);
end
