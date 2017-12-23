function [Xh] = applyBC(Xh, BCNumber, n)
% applies boundary conditions to vector by extending it
% data:
%     Xh          vector containing the solution at current time/ RK step
%     BCNumber    string specifying which boundary condition to apply 
%                 ("Periodic" or "Open")
% output:
%     Xh    vecor of current solution (by columns) extended at its beginning and
%           end with the selected BC
    
    N = size(Xh,1);
    switch BCNumber
        case 0
            Xh = [Xh(N-n+1:N,:); Xh; Xh(1:n,:)];
        case 1
            Xh = [repmat(Xh(1,:),n,1); Xh; repmat(Xh(end,:),n,1)];
    end
end
