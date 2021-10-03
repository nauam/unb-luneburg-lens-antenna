function diagram(diag,fig_num,func,fig_title,windows,xstep,ystep)
double precision; format long; 
    figure(fig_num);
    if diag == "polar"
        polarplot(0:pi/500:2*pi,func,'Color',[0 0 0],'Linewidth',2);
        pax = gca;
        pax.ThetaLim = [-180 180];
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
    else %diag == "cart"
        plot((0:pi/500:2*pi)*180/pi,20*log10(func),'Color',[0 0 0],'Linewidth',2);
        xlabel('Theta [deg]'); 
        ylabel('dB20(rETotal)');
        xmax = 360;                             xmin = 0;              
        ymax = ceil(max(20*log10(func))/10)*10; ymin = ymax - windows; 
        axis([xmin xmax ymin ymax]); 
        grid on;
        xticks(xmin:xstep:xmax);
        yticks(ymin:ystep:ymax);
    end
    title(fig_title);
end