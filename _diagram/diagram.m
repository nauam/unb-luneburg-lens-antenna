function diagram(diag,num,func)
double precision; format long; 
    figure(num);
    title('Far field pattern');
    if diag == "polar"
        polarplot(0:pi/500:2*pi,func,'Color',[0 0 0],'Linewidth',2);
        pax = gca;
        pax.ThetaLim = [-180 180];
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
    else %diag == "cart"
        plot((0:pi/500:2*pi)*180/pi,20*log10(func),'Color',[0 0 0],'Linewidth',2);
        xlabel('Theta [deg]');
        ylabel('dB20normalize(rETotal)');
        axis([0 360 -40 0]); 
        grid on;
    end
end