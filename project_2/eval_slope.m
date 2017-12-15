function [ slope ] = eval_slope( ul, ur, SlopeNumber)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    N = length(ul);
    switch SlopeNumber
        case 0
            slope = zeros(N,1);
        case 1
            slope = minmod([ul,ur]);
        case 2
            slope = minmod([(ul+ur)/2 2*ul 2*ur]);
    end

end

