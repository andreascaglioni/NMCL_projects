function [ERR] = computeL1Error(hExa, mExa, hh, mh, xx, T)
    N = length(xx)-1;
    err = zeros(1,N);
    for j=1:N
        hExaFinal =@(x) abs(hExa(x,T)-hh(j));
        mExaFinal =@(x) abs(mExa(x,T)-mh(j));
        err(j) = integral(hExaFinal,xx(j),xx(j+1),'absTol', 1e-14)...
               + integral(mExaFinal,xx(j),xx(j+1),'absTol', 1e-14);
    end
    ERR = sum(err);
end

