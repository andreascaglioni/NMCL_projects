    function [xx, Xh, uh, ERR]=main(dx)
% the "main" function fo the assignment which executes the task and calls
% the other functions. The only data is the meshs spacing, however all
% other parametnes can be chosen in the section "choose" (boundary 
% conditions, initial conditions, source term, exact solution, numerical 
% flux) (see comments) and "domain and discretization" (CFL number, final
% time).
% data:
%     dx      spacing for space discretization
% returns:    
%     xx      vector containig centerpoints of cells
%     Xh      vector containing the solution at final time T (by columns)
%     uh      numerical solution for speed (ratio discharge/height) at 
%             final time
%     ERR     L1 error of the numerical solution  at final time (n.b. 
%     relevant only if the analytical solution is available
%% Choose
BCNumber = 1;   % 0:periodic;   1:open
ICNumber = 3;   % 0: ex. 1.1;   1,2: ex. 1.2(a) and 1.2(b);     3: ex. 1.4
SourceNumber = 1;  % 0: ex: 1.1;   1: no source
ExaNumber = 1;  % 0:ex. 1.1;   1: no exact sol available (put it to 0)
FluxNumber = 0; % 0:LF;     1:Roe
SlopeNumber = 1;    %0: zero;   1: MinMod;  2:MUSCL
%% physical parameters
g = 1;
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
        m0 = @(x) (-1.5*(x<1.));
end
%% exact solution
switch ExaNumber
    case 0
        hExa = @(x,t) h0(x-t);
        mExa = @(x,t) u * hExa(x,t);
    case 1
        hExa = @(x,t) 0*x+0*t;
        mExa = @(x,t) 0*x+0*t;
end
%% domain and discretization
cc = 0:dx:2;    % cells' boundaries
xx = dx/2:dx:2-dx/2;    % cells' centerpoints
N = length(xx); % number of cells
CFL = 0.5;
T = 0.5;
%% discrete IC (integrated)
hh = zeros(N, 1);
mh = zeros(N, 1);
for j = 1:N
    hh(j) = integral(h0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
    mh(j) = integral(m0,cc(j),cc(j+1),'AbsTol',1e-14)/dx;
end
%% velocity and solution vectors
uh = mh./hh;
Xh = [hh, mh];
%% solve
time = 0.;
nIter = 0;
while time < T
    %compute new timestep
    maxVel = max(abs(uh)+(g*Xh(:,1)).^0.5);
    k = CFL*dx/maxVel;
    if(time+k>T)
        k = T-time;
    end
    % Update solution
    rhs = EulerSLrhs1D(Xh, xx, dx, time, maxVel, BCNumber, FluxNumber, SlopeNumber, SourceNumber,cc); 
    X1 = Xh + k*rhs;
    rhs = EulerSLrhs1D(X1, xx, dx, time+k, maxVel, BCNumber, FluxNumber, SlopeNumber, SourceNumber,cc); 
    X2 = (3*Xh + X1 + k*rhs)/4;
	rhs = EulerSLrhs1D(X2, xx, dx, time+0.5*k, maxVel, BCNumber, FluxNumber, SlopeNumber, SourceNumber,cc); 
    Xh = (Xh + 2*X2 + 2*k*rhs)/3;
    uh = Xh(:,2)./Xh(:,1);
    %plot sol at every timestep
    if(mod(nIter,1)==0)
      plotSolAtTime(xx, Xh,  @(x)hExa(x,time), @(x)mExa(x,time));
%       plotErrAtTime(@(x)hExa(x,time), @(x)mExa(x,time), Xh, xx, dx, cc, time);
%       press = waitforbuttonpress;
      pause(0.0001);
    end
   time = time+k;
   nIter = nIter +1;
end
%% plot final solution
figure
plotSolAtTime(xx, Xh,  @(x)hExa(x,T), @(x)mExa(x,T));
%% compute error at final time (in L^1)
ERR = computeL1Error(@(x)hExa(x,time), @(x)mExa(x,time), Xh, cc, dx);
end
