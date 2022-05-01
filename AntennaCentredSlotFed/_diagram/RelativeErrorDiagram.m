function RelativeErrorDiagram(n, nfigure, antenna, settings, lgnd)
    addpath(genpath('_function')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    double precision; format long;

    for settings_ = settings
        [Nlens, k0r, omega, r, ang, epsilonR, source, ~] = Structure(antenna, settings_);
        [N, R, Theta, Epsilon, Mu] = LuneburgLensStructure(Nlens, k0r, omega, r, ang, epsilonR);

        f = waitbar(0, strcat('n = ', num2str(n)));
        for n_ = 1:n + 1
            [X(1:n_, 1:N, n_), ~] = RegularizationFunction(n_, N, omega, source, R, Theta, Epsilon, Mu);
            waitbar(n_ / n, f, sprintf(strcat('n_ = ', num2str(n_), '` de `', num2str(n + 1), '` e `', lgnd, '` = `', num2str(Nlens), '`')))
        end
        delete(f)
        
        index = find(settings == settings_);
        for n_ = 1:n
            erro(n_, index) = norm(X(1:n_, :, n_)' - X(1:n_, :, n_ + 1)', Inf) / norm(X(1:n_, :, n_ + 1)', Inf);
        end
        legend_str{index} = strcat(lgnd, '=', num2str(Nlens));

    end

    %diagram
    figure(nfigure)
    semilogy(1:n, erro, 'Linewidth', 2);
    hold on
    grid on;
    legend(legend_str, 'Location', 'northwest', 'NumColumns', index);
    xlabel('Tuncamento n');
    ylabel('Erro Relativo e_r(n)');

end
