function [hh, mh] = applyBC(hh, mh, bc)
% applies boundary conditions to vector by extending it
% data:
%     hh    vecor of first unknown (height)
%     mh    vecor of second unknown (discharge)
%     bc    string specifying which boundary condition to apply ("Periodic"
%           or "Open")
% output:
%     hh    vecor of first unknown (height) extended at its beginning and
%           end with the selected BC
%     mh    vecor of second unknown (discharge) extended at its beginning and
%           end with the selected BC
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
