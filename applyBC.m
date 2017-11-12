function [hh, mh] = applyBC(Yh, bc)
    N = length(Yh)/2;
    hh = Yh(1:N,1); 
    mh = Yh(N+1:N*2,1);
    hh1 = hh(1);
    mh1 = mh(1);
    hhN = hh(end); 
    mhN = mh(end);
    switch bc
        case 'Periodic'
            hh = [hhN; hh; hh1];
            mh = [mhN; mh; mh1];
        case 'Open'
            hh = [hh1; hh; hhN];
            mh = [mh1; mh; mhN];
    end
end
