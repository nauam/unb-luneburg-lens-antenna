function S = SFunction(u_, n_, Theta)
    double precision;
    format long;

    S = zeros(length(Theta));

    for i = 1:length(Theta)
        theta = Theta(i, i);

        if theta == 0
            S(i, i) = 0;
        else
            S(i, i) = 2 / pi * (QFunction(u_, n_, theta) - QFunction(u_, 0, theta) * QFunction(0, n_, theta) / QFunction(0, 0, theta));
        end

    end

end
