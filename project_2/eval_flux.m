function [ Fh ] = eval_flux(Xhl, Xhr, FluxNumber, f1, f2, g,maxVel)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    switch FluxNumber
        case 0
            Fh = LFNumericalFlux(Xhl, Xhr, f1, f2, g, maxVel);
        case 1
            Fh = RoeNumericalFlux(Xhl, Xhr, f1, f2, g, maxVel);
end

