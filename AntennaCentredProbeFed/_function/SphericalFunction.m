function [Z1, Z3, K1, K3] = SphericalFunction(Z, n)
    double precision;
    format long;

    N = length(Z);
    Z1 = zeros(N, N, n);
    Z3 = zeros(N, N, n);
    K1 = zeros(N, N, n);
    K3 = zeros(N, N, n);

    n_ = (1:n);
    n_half = (1:(n + 1)) + 1/2;

    for i = 1:N
        z = Z(i, i);

        jn = besselj(n_half, z) .* sqrt(pi / (2 * z));
        yn = bessely(n_half, z) .* sqrt(pi / (2 * z));
        hn = jn + 1i * yn;

        Z1(i, i, :) = jn(n_);
        Z3(i, i, :) = hn(n_);
        K1(i, i, :) = jn(n_) .* (n_ + 1) / z - jn(n_ + 1);
        K3(i, i, :) = hn(n_) .* (n_ + 1) / z - hn(n_ + 1);
    end

end
