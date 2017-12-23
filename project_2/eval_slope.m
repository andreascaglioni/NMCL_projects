function [ slope ] = eval_slope( ul, ur, SlopeNumber)
% allows to choose what slope to use and computes it
% Data:
%    ul, ur      vectors of vaues to cmpute slope
%    SlopeNumber number to choose which slope to compute
%   Detailed explanation goes here
    N = length(ul);
    switch SlopeNumber
        case 0
            slope = zeros(N,2);
        case 1
            slope(:,1) = minmod([ul(:,1), ur(:,1)]);
            slope(:,2) = minmod([ul(:,2), ur(:,2)]);
        case 2
            slope(:,1) = minmod([(ul(:,1)+ur(:,1))/2. , 2.*ul(:,1), 2.*ur(:,1)]);
            slope(:,2) = minmod([(ul(:,2)+ur(:,2))/2. , 2.*ul(:,2), 2.*ur(:,2)]);
    end
end

