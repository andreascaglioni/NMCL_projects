function [] = plotSolAtTime(cc, Xh, hExa, mExa)
% plots given functions of 1 varaible that take values in R^2 compunentwise
% data: 
%     cc    points for x component of plot
%     Xh    vector of values of function to plot (compunents by columns)
%     hExa  function handle hExa(x) depending on space for exact
%           values of first compunent (at given time)
%     mExa  function handle mExa(x,t) depending on space and time for exact
%           values of second component (at given time)
    hh = Xh(:,1);
    mh = Xh(:,2);
    
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
end

