function []=plotErrAtTime(hExa, mExa, Xh, xx, dx, cc, time)
    N = length(cc)-1;
    err = zeros(2,N);
    for j=1:N
        % error using exact solution integrated with high-order scheme
%         hDiff =@(x) abs(hExa(x)-Xh(j,1));
%         mDiff =@(x) abs(mExa(x)-Xh(j,2));
%         err(j) = integral(mDiff,cc(j),cc(j+1),'absTol', 1e-14) +...
%                  integral(mDiff,cc(j),cc(j+1),'absTol', 1e-14);

        %using cell averages of exact solution
        hExaAver(j) = integral(hExa, cc(j), cc(j+1), 'absTol', 1e-14)/dx;
        mExaAver(j) = integral(mExa, cc(j), cc(j+1), 'absTol', 1e-14)/dx;
        err(1,j) = dx*abs(hExaAver(j)-Xh(j,1));
        err(2,j) = dx* abs(mExaAver(j)-Xh(j,2));
    end
%     figure
    subplot(2,1,1)
    plot(cc(1:end-1), err(1,:));
    subplot(2,1,2);
    plot(cc(1:end-1), err(2,:));
    title(num2str(time));
    ERR = sum(err(:));%sqrt(sum(err.^2));
    
end