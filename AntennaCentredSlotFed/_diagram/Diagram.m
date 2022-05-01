function Diagram(nFigure, dTh, titleFigure, Axis, color)
    double precision;
    format long;

    figure(nFigure);

    %Far Field Diagram
    if Axis(1) == "dB"
        CPTh = 20 * log10(abs(dTh) / max(abs(dTh)));
    else % AxisF(1) == "norm"
        CPTh = abs(dTh) / max(abs(dTh));
    end

    mn = str2double(Axis(3));
    mx = str2double(Axis(4));
    step = str2double(Axis(5));

    CPTh(251) = nan;
    CPTh(751) = nan;

    if Axis(2) == "polar"
        polarplot(0:pi / 500:2 * pi, CPTh, 'Color', color, 'Linewidth', 2);
        rlim([mn mx])
        pax = gca;
        pax.ThetaLim = [-180 180];
        pax.ThetaDir = 'clockwise';
        pax.ThetaZeroLocation = 'top';
    else %diagram == "cart"
        plot((0:pi / 500:2 * pi) * 180 / pi, CPTh, 'Color', color, 'Linewidth', 2);
        xlabel('Theta [deg]');
        ylabel('dB20(rETotal)');
        grid on;

        xmax = 360;
        xmin = 0;
        xstep = 60;

        xticks(xmin:xstep:xmax);
        yticks(mn:step:mx);
        axis([xmin xmax mn mx]);
    end

    title(titleFigure);
end
