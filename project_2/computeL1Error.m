function [ERR] = computeL1Error(hExa, mExa, Xh, xx, dx)
% computes L1 norm of the difference of numerical solution and given exact solution
% data:
%      hExa  function handle hExa(x) depending on space for exact
%            values of hh (at given time)
%      mExa  function handle mExa(x,t) depending on space and time for exact
%            values of mh (at given time)
%      hh    vector of values of numerical solution (height)
%      mh    vector of values of numerical solution (discharge)
%      xx    vactor of x values between which to integrate
% returns:
%     ERR     real number (the 1 error)
    N = length(xx)-1;
    hh = Xh(:,1);
    mh = Xh(:,2);
    err = zeros(1,N);
    for j=1:N
        % error using exact solution integrated with high-order scheme
        hExaFinal =@(x) abs(hExa(x)-hh(j));
        mExaFinal =@(x) abs(mExa(x)-mh(j));
        err(j) = (integral(hExaFinal,xx(j),xx(j+1),'absTol', 1e-14)+...
                  integral(mExaFinal,xx(j),xx(j+1),'absTol', 1e-14));
        
        % using cell averages of exact solution
%         hExaAver = integral(hExa,xx(j),xx(j+1),'absTol', 1e-14)/dx;
%         mExaAver = integral(mExa,xx(j),xx(j+1),'absTol', 1e-14)/dx;
%         err(j) = dx*(abs(hExaAver - hh(j)) + abs(mExaAver - mh(j)));
    end
    ERR = sum(err);
end

