function [] = plotSOlAtTIme(cc, hh, hExa, mh, mExa)
% plots given functions of 1 varaible that take values in R
% data: 
%     cc    points for x component of plot
%     hh    vector of values of a function (height) to plot
%     hExa  function handle hExa(x) depending on space for exact
%           values of hh (at given time)
%     mh    vector of values of a function (discharge) to plot
%     mExa  function handle mExa(x,t) depending on space and time for exact
%           values of mh (at given time)
    uh = mh./hh;
    subplot(3,1,1)
    plot(cc, hh, cc ,hExa(cc));
    grid on; title 'height';
    subplot(3,1,2)
    plot(cc, mh, cc ,mExa(cc));
    grid on; title 'discharge';
    subplot(3,1,3)
    plot(cc, uh); grid on; title 'speed';
end

