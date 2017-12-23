function [ERR] = computeL1Error(hExa, mExa, Xh, cc, dx)
% computes L1 norm of the difference of numerical solution and given exact solution
% data:
%     hExa   function handle hExa(x) depending on space for first component
%            of exact solution
%     mExa   function handle hExa(x) depending on space for second component
%            of exact solution
%     Xh     vector containing the solution at final time T (by columns)
%     cc     vector of x values between which to integrate
%     dx     mesh spacing
% returns:
%     ERR     real number (the L^1 error)
    N = length(cc)-1;
    err = zeros(N,2);
    XExaAver = zeros(N,2);
    for j=1:N
        %using cell averages of exact solution
        XExaAver(j,1) = integral(hExa, cc(j), cc(j+1), 'absTol', 1e-14)/dx;
        XExaAver(j,2) = integral(mExa, cc(j), cc(j+1), 'absTol', 1e-14)/dx;
        err(j,1) = dx * abs(XExaAver(j,1)-Xh(j,1));
        err(j,2) = dx * abs(XExaAver(j,2)-Xh(j,2));
    end
    ERR = sum(err(:));
end

