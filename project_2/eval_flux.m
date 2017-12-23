function [ Fh ] = eval_flux(u, v, FluxNumber,maxVel)
% allows to choose which slope to use an computes it
% Data:
%     u,v         vectors containig the values to use to compute the slope
%     FluxNumber  number to choose the numerical flux
%     maxVel      maximum velocity of the system
% Returns:
%     Fh          vaactor containing (by columns) values of flux
%% choose and compute flux
    switch FluxNumber
        case 0
            Fh = LFNumericalFlux(u, v, maxVel);
        case 1
            Fh = RoeNumericalFlux(u, v, maxVel);
    end
end

