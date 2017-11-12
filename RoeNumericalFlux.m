function [F1Ext, F2Ext] = LFNumericalFlux(hh, mh, f1, f2, g);
    N = length(hh)-2;
    z1 = 0.5*(hh(2:end).^0.5 + hh(1:end-1).^0.5);
    z1s = 0.5*(hh(2:end) + hh(1:end-1));
    z2 = 0.5*(mh(2:end)./((hh(2:end).^0.5)) + mh(1:end-1)./((hh(1:end-1).^0.5)));
    w = z2./z1;
    F1Ext = zeros(N,1);
    F2Ext = zeros(N,1);
    for i = 1:N+1
        lambdaMax = w(i)+(g*z1s(i))^0.5;
        lambdaMin = w(i)-(g*z1s(i))^0.5;
        vMax = 1/(1+(w(i)+(g*z1s(i))^0.5)^2)^0.5*[1; w(i)+(g*z1s(i))^0.5];
        vMin = 1/(1+(-w(i)+(g*z1s(i))^0.5)^2)^0.5*[1;w(i)-(g*z1s(i))^0.5];
        SRoe = [vMin, vMax];
        LRoe = [abs(lambdaMin), 0; 0, abs(lambdaMax)];
        ARoe = SRoe\LRoe*SRoe;        
        % numerical flux
        F1Ext(i) = 0.5*(f1(hh(i+1),mh(i+1)) + f1(hh(i),mh(i))) - 0.5*(ARoe(1,:)*[hh(i+1)-hh(i); mh(i+1)-mh(i)]);
        F2Ext(i) = 0.5*(f2(hh(i+1),mh(i+1)) + f2(hh(i),mh(i))) - 0.5*(ARoe(2,:)*[hh(i+1)-hh(i); mh(i+1)-mh(i)]);
    end
end

