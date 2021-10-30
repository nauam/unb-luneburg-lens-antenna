function diagram(diag,fig_num,func,fig_title,var1,var2)
double precision; format long; 
    figure(fig_num);
    if diag == "polar"
        polarplot(0:pi/500:2*pi,func,'Color',[0 0 0],'Linewidth',2);
        rlim([var1 var2])
        pax = gca;
        pax.ThetaLim = [-180 180];
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'right';
    else %diag == "cart"
        plot((0:pi/500:2*pi)*180/pi,func,'Color',[0 0 0],'Linewidth',2);
        xlabel('Theta [deg]'); 
        ylabel('dB20(rETotal)');
        xmax = 360;                 xmin = 0;                       xstep = ceil((xmax-xmin)/var1);
        ymax = ceil(max(func)/5)*5; ymin = ceil(min(func)/5)*5; ystep = ceil((ymax-ymin)/var2);
        axis([xmin xmax ymin ymax]); 
        grid on;
        xticks(xmin:xstep:xmax);
        yticks(ymin:ystep:ymax);
    end
    title(fig_title);
end