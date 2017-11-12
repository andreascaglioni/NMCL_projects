function [] = plotSOlFInalTIme(cc, hh, hExa, mh, mExa, uh, time)
    %figure
    subplot(3,1,1)
    plot(cc, hh, cc ,hExa(cc,time));
    grid on; title 'height';
    subplot(3,1,2)
    plot(cc, mh, cc ,mExa(cc,time));
    grid on; title 'discharge';
    subplot(3,1,3)
    plot(cc, uh); grid on; title 'speed';
end

