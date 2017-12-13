function [Uh] = applyBC(Uh, bcNumber)
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
%     hh1 = hh(1);
%     mh1 = mh(1);
%     hhN = hh(end); 
%     mhN = mh(end);
    Uh1 = Uh(1,:);
    N = size(Uh,1);
    UhN = Uh(N,:);
    switch bcNumber
        case 0
              Uh = [UhN; Uh; Uh1];
        case 1
              Uh = [Uh1; Uh; UhN];
    end
end
