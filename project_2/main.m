    function [xx, Xh, uh, ERR]=main(dx)
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
FluxNumber = 0; % 0:LF;     1:Roe
SlopeNumber = 0;    %0: zero;   1: MinMod;  2:MUSCL
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
xx = 0:dx:2;    %cells boundaries
cc = dx/2:dx:2-dx/2;    %cells' centerpoints
N = length(cc); % number of cells
CFL = 0.5;
T = 2.;
%% discrete IC (integrated)
hh0 = zeros(N, 1);      %stores discrete solution in a 1 dim vector (for every t)
mh0 = zeros(N, 1);
for j = 1:N
    hh0(j) = integral(h0,xx(j),xx(j+1),'AbsTol',1e-14)/dx;
    mh0(j) = integral(m0,xx(j),xx(j+1),'AbsTol',1e-14)/dx;
end
%% solution vectors
hh = hh0;  %initialize solution to IC
mh = mh0;
uh = mh0./hh0;
Xh = [hh, mh];
% plotSolAtTime(xx, hh, @(x)hExa(x,0.), mh, @(x)mExa(x,0.));
% press  = waitforbuttonpress;
%% solve
time = 0.;
% figure;
while time < T
    %compute new timestep
    maxVel = max(abs(uh)+(g*Xh(:,1)).^0.5);
    k = CFL*dx/maxVel;
    if(time+k>T)
        time = T-k;
    end
    
    % Update solution
    rhs  = eval_rhs(Xh, BCNumber, FluxNumber, SlopeNumber, f1, f2, g, S1Prim, S2Prim, xx, dx, time, maxVel); 
    X1 = Xh + k*rhs;
    rhs  = eval_rhs(X1, BCNumber, FluxNumber, SlopeNumber, f1, f2, g, S1Prim, S2Prim, xx, dx, time, maxVel); 
    X2 = (3*Xh + X1 + k*rhs)/4;
	rhs  = eval_rhs(X2, BCNumber, FluxNumber, SlopeNumber, f1, f2, g, S1Prim, S2Prim, xx, dx, time, maxVel); 
    Xh  = (Xh + 2*X2 + 2*k*rhs)/3;

    uh = Xh(:,2)./Xh(:,1);
    
    time = time+k;

    %plot sol at every timestep
%    plotSolAtTime(xx, hh, @(x)hExa(x,time), mh, @(x)mExa(x,time));
%    press = waitforbuttonpress;
%   pause(0.001);
end
%% plot final solution
figure
plotSolAtTime(cc, Xh,  @(x)hExa(x,T), @(x)mExa(x,T));
%% compute error at final time (in L^1)
ERR = computeL1Error(hExa, mExa, Xh, cc, T);
end
