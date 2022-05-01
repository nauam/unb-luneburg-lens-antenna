function [dTh, dPh] = FarFieldFunction(m, n, a3e0mn, a3o0mn, b3e0mn, b3o0mn, phi)
    double precision;
    format long;

    t1mn = zeros(1001, n);
    t2mn = zeros(1001, n);
    cos_t = cos(0:pi / 500:2 * pi);
    sin_t = sin(0:pi / 500:2 * pi);

    for n_ = 1:n
        Pmn = legendre(n_, cos_t);
        dPmn = LegendreDerivativeFunction(n_, cos_t);
        t1mn(:, n_) = m * Pmn(m + 1, :) ./ sin_t;
        t2mn(:, n_) = dPmn(m + 1, :) .* sin_t;
    end

    mip1 = diag((-1i).^((1:n) + 1));

    dTh = (t1mn * mip1 * a3o0mn - 1i * t2mn * mip1 * b3o0mn) * cos(m * phi) - (t1mn * mip1 * a3e0mn + 1i * t2mn * mip1 * b3e0mn) * sin(m * phi);
    dPh = (t2mn * mip1 * a3e0mn + 1i * t1mn * mip1 * b3e0mn) * cos(m * phi) + (t2mn * mip1 * a3o0mn - 1i * t1mn * mip1 * b3o0mn) * sin(m * phi);
end
