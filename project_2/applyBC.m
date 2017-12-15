function [Xh] = applyBC(Xh, BCNumber, n)
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
    
    N = size(Xh,1);
    Xbeg = Xh(1:n,:);
    Xend = Xh(N-n+1:N,:);
    switch BCNumber
        case 0
            Xh = [Xend; Xh; Xbeg];
        case 1
            Xh = [Xbeg; Xh; Xend];
    end
end
