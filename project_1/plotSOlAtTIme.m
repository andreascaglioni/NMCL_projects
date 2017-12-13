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
    
    x0=10;
    y0=10;
    width=550;
    height=700;
    
    subplot(2,1,1)
    plot(cc, hh, cc ,hExa(cc));
    set(gcf,'units','points','position',[x0,y0,width,height]);
    grid on;
    title 'Height';
    set(gca,'FontSize',16);
    subplot(2,1,2)
    plot(cc, mh, cc ,mExa(cc));
    set(gcf,'units','points','position',[x0,y0,width,height]);
    grid on;
    title 'Discharge';
    set(gca,'FontSize',16);
%     subplot(3,1,3)
%     plot(cc, uh);
%     set(gcf,'units','points','position',[x0,y0,width,height]);
%     grid on;
%     title 'Speed';
%     set(gca,'FontSize',16);
end

