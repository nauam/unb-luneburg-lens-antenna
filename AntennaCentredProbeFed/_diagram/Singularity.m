function Singularity(nFigure, titleFigure, Gbn, n)
    double precision;
    format long;

    figure(nFigure);
    single = zeros(1, n - 1);

    for n_ = 1:n - 1
        single(n_) = norm(diag(Gbn(:, :, n_)) - diag(Gbn(:, :, n_ + 1)), Inf) / norm(diag(Gbn(:, :, n_ + 1)), Inf);
    end

    %diagram
    semilogy(1:n - 1, single, 'Linewidth', 2);
    hold on
    grid on;
    xlabel('Tuncamento n');
    ylabel('Erro Relativo e_r(n)');
    title(titleFigure);
end
