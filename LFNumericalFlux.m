function [F1Ext, F2Ext] = LFNumericalFlux(hh, mh, f1, f2, g)
    uh = mh./hh;
    lMax =max((uh(1:end)+(g*hh(1:end)).^0.5));	% maximum wave speed (vector)
    F1Ext = 0.5*(f1(hh(2:end),mh(2:end)) + f1(hh(1:end-1),mh(1:end-1)))-0.5*lMax*(hh(2:end) - hh(1:end-1));
    F2Ext = 0.5*(f2(hh(2:end),mh(2:end)) + f2(hh(1:end-1),mh(1:end-1)))-0.5*lMax*(mh(2:end) - mh(1:end-1));
end

