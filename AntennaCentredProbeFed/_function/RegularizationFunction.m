function [XuN, b3nN, Gbn] = RegularizationFunction(n, N, omega, source, R, Theta, Epsilon, Mu)
    double precision;
    format long;

    b3nN = zeros(n, N);
    Ykbn_ = zeros(N, n);
    Yzbn_ = zeros(N, n);
    Gbn = zeros(N, N, n);
    CZbn = zeros(N, N, n);
    Sun_ = zeros(N, n, N, n);
    SunGbn_ = zeros(N, n, N, n);

    I = eye(N * n);
    I_N = eye(N);
    Is = diag(ones(1, (N - 1)), 1);
    Ii = diag(ones(1, (N - 1)), -1);
    i0 = 1:N;
    i1 = 2:N + 1;
    d0 = double(source(1) == i0 - 1)';
    d1 = double(source(1) == i1 - 1)';

    Eti0 = sqrt(Mu(i0, i0) / Epsilon(i0, i0));
    Eti1 = sqrt(Mu(i1, i1) / Epsilon(i1, i1));
    Ki0 = omega * sqrt(Mu(i0, i0) * Epsilon(i0, i0));
    Ki1 = omega * sqrt(Mu(i1, i1) * Epsilon(i1, i1));
    ks = Ki1(source(1), source(1));

    fk = ones(1, n);
    fz = (2 * (1:n) + 1) ./ ((1:n) .* ((1:n) + 1));
    OmgB =- inv(R * (Ki0 / Eti0 + Ki1 / Eti1) / 2);

    [Z1i0, Z3i0, K1i0, K3i0] = SphericalFunction(Ki0 * R, n);
    [Z1i1, Z3i1, K1i1, K3i1] = SphericalFunction(Ki1 * R, n);
    [Z1s, Z3s, ~, ~] = SphericalFunction(ks * R(source(2), source(2)), n);

    for n_ = 1:n
        B1s = Z1s(1, 1, n_) * (2 * n_ + 1) * d1;
        B3s = Z3s(1, 1, n_) * (2 * n_ + 1) * d0;

        K1Dbn = K1i1(:, :, n_);
        K1Ebn = K1i0(:, :, n_) * Ii;
        K1Fbn = K1i0(:, :, n_);
        K1Tbn = K1Dbn - K1Ebn;

        K3Dbn = K3i0(:, :, n_);
        K3Ebn = K3i1(:, :, n_) * Is;
        K3Fbn = K3i1(:, :, n_);
        K3Tbn = K3Dbn - K3Ebn;

        Z1Dbn = Z1i1(:, :, n_) / Eti1;
        Z1Ebn = Z1i0(:, :, n_) / Eti0 * Ii;
        Z1Fbn = Z1i0(:, :, n_) / Eti0;
        Z1Tbn = Z1Dbn - Z1Ebn;

        Z3Dbn = Z3i0(:, :, n_) / Eti0;
        Z3Ebn = Z3i1(:, :, n_) / Eti1 * Is;
        Z3Fbn = Z3i1(:, :, n_) / Eti1;
        Z3Tbn = Z3Dbn - Z3Ebn;

        CKbn_ = K3Dbn + K1Ebn / K1Tbn * K3Tbn;
        CZbn_ = Z3Tbn - Z1Tbn / K1Tbn * K3Tbn;
        Ckbn_ = K1Ebn / K1Tbn * K3Fbn * B1s - (eye(N) + K1Ebn / K1Tbn) * K1Fbn * B3s;
        Czbn_ = (Z3Fbn - Z1Tbn / K1Tbn * K3Fbn) * B1s - (Z1Fbn - Z1Tbn / K1Tbn * K1Fbn) * B3s;

        Yzbn_ (:, n_) = Czbn_ / fz(n_);
        Ykbn_ (:, n_) = (OmgB \ Ckbn_) / fk(n_);
        Gbn(:, :, n_) = I_N - fz(n_) / fk(n_) * (OmgB \ CKbn_) / CZbn_;
        CZbn(:, :, n_) = CZbn_;
    end

    for u_ = 1:n

        for n_ = 1:n
            Sun_ (:, u_, :, n_) = SFunction(u_, n_, Theta);
            SunGbn_(:, u_, :, n_) = SFunction(u_, n_, Theta) * Gbn(:, :, n_);
        end

    end

    Ykbn = reshape(Ykbn_, [N * n 1]);
    Yzbn = reshape(Yzbn_, [N * n 1]);
    Sun = reshape(Sun_, [N * n N * n]);
    SunGbn = reshape(SunGbn_, [N * n N * n]);

    Ybn = Yzbn + Sun * (Ykbn - Yzbn);
    Xu_ = (I - SunGbn) \ Ybn;

    XuN = reshape(Xu_, [N n])';

    for n_ = 1:n
        b3nN(n_, :) = fz(n_) * (CZbn(:, :, n_) \ XuN(n_, :)');
    end

end
