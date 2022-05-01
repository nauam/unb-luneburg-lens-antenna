function RadiatedPowerDiagram(n, nfigure, antenna, settings, range, step)
    addpath(genpath('_function')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    double precision; format long;

    charge = range/step;
    
    for settings_ = settings
        f = waitbar(0, strcat('charge = ', num2str(charge)));
        for charge_ = 1:(charge + 1)
            k0r = step * charge_;
            [Nlens, ~, omega, r, ang, epsilonR, source, ~] = Structure(antenna, settings_);
            [N, R, Theta, Epsilon, Mu] = LuneburgLensStructure(Nlens, k0r, omega, r, ang, epsilonR);
            [~, b3nN, ~] = RegularizationFunction(n, N, omega, source, R, Theta, Epsilon, Mu);
            Pr_norm(charge_, settings_) = RadiatedPowerFunction(n, 0, 0, 0, b3nN(:, 1));
            waitbar(charge_ / charge, f, sprintf(strcat('charge_ = ', num2str(charge_), '` de `', num2str(charge + 1), '` e N = `', num2str(Nlens), '`')))
        end
        delete(f)
        indice = find(Pr_norm(1:charge + 1, settings_) == max(Pr_norm(1:charge + 1, settings_))) * .05;
        legend_str{settings_} = strcat('N=', num2str(Nlens), '` k0r=', num2str(indice));
    end

    %diagram
    figure(nfigure)
    semilogy((0:charge) * step, abs(Pr_norm), 'Linewidth', 2);
    hold on
    grid on;
    legend(legend_str, 'Location', 'northwest', 'NumColumns', 5);
    ylabel('Potência irradiada normalizada');
    xlabel('Frequência normalizada');

end
