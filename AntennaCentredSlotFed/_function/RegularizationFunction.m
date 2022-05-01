function [XauN, XbuN, a3nN, b3nN, Gan, Gbn] = RegularizationFunction(n, N, omega, source, R, Theta, Epsilon, Mu)
    double precision;
    format long;

    a3nN = zeros(n, N);
    Yzan_ = zeros(N, n);
    Ykan_ = zeros(N, n);
    Gan = zeros(N, N, n);
    CKan = zeros(N, N, n);

    b3nN = zeros(n, N);
    Ykbn_ = zeros(N, n);
    Yzbn_ = zeros(N, n);
    Gbn = zeros(N, N, n);
    CZbn = zeros(N, N, n);

    Sun_ = zeros(2 * N, n, 2 * N, n);
    SunGabn_ = zeros(2 * N, n, 2 * N, n);

    I = eye(2 * N * n);
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

    fk = (2 * (1:n) + 1) ./ ((1:n) .* ((1:n) + 1));
    fz = 1 ./ ((1:n) .* ((1:n) + 1));
    OmgA =- inv(R \ (inv(Ki0 * Eti0) + inv(Ki1 * Eti1)) / 2);
    OmgB =- inv(R * (Ki0 / Eti0 + Ki1 / Eti1) * 2);

    [Z1i0, Z3i0, K1i0, K3i0] = SphericalFunction(Ki0 * R, n);
    [Z1i1, Z3i1, K1i1, K3i1] = SphericalFunction(Ki1 * R, n);
    [Z1s, Z3s, K1s, K3s] = SphericalFunction(ks * R(source(2), source(2)), n);

    for n_ = 1:n
        A1s = K1s(1, 1, n_) * (2 * n_ + 1) / (n_ * (n_ + 1)) * d1;
        A3s = K3s(1, 1, n_) * (2 * n_ + 1) / (n_ * (n_ + 1)) * d0;
        B1s = Z1s(1, 1, n_) * (2 * n_ + 1) / (n_ * (n_ + 1)) * d1;
        B3s = Z3s(1, 1, n_) * (2 * n_ + 1) / (n_ * (n_ + 1)) * d0;

        Z1Dan = Z1i1(:, :, n_);
        Z1Ean = Z1i0(:, :, n_) * Ii;
        Z1Fan = Z1i0(:, :, n_);
        Z1Tan = Z1Dan - Z1Ean;

        Z3Dan = Z3i0(:, :, n_);
        Z3Ean = Z3i1(:, :, n_) * Is;
        Z3Fan = Z3i1(:, :, n_);
        Z3Tan = Z3Dan - Z3Ean;

        K1Dan = K1i1(:, :, n_) / Eti1;
        K1Ean = K1i0(:, :, n_) / Eti0 * Ii;
        K1Fan = K1i0(:, :, n_) / Eti0;
        K1Tan = K1Dan - K1Ean;

        K3Dan = K3i0(:, :, n_) / Eti0;
        K3Ean = K3i1(:, :, n_) / Eti1 * Is;
        K3Fan = K3i1(:, :, n_) / Eti1;
        K3Tan = K3Dan - K3Ean;

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

        CZan_ = Z3Dan + Z1Ean / Z1Tan * Z3Tan;
        CKan_ = K3Tan - K1Tan / Z1Tan * Z3Tan;
        CKbn_ = K3Dbn + K1Ebn / K1Tbn * K3Tbn;
        CZbn_ = Z3Tbn - Z1Tbn / K1Tbn * K3Tbn;

        Czan_ = Z1Ean / Z1Tan * Z3Fan * A1s - (eye(N) + Z1Ean / Z1Tan) * Z1Fan * A3s;
        Ckbn_ = K1Ebn / K1Tbn * K3Fbn * B1s - (eye(N) + K1Ebn / K1Tbn) * K1Fbn * B3s;
        Ckan_ = (K3Fan - K1Tan / Z1Tan * Z3Fan) * A1s - (K1Fan - K1Tan / Z1Tan * Z1Fan) * A3s;
        Czbn_ = (Z3Fbn - Z1Tbn / K1Tbn * K3Fbn) * B1s - (Z1Fbn - Z1Tbn / K1Tbn * K1Fbn) * B3s;

        Ykan_ (:, n_) = Ckan_ / fk(n_);
        Yzbn_ (:, n_) = Czbn_ / fz(n_);
        Yzan_ (:, n_) = (OmgA \ Czan_) / fz(n_);
        Ykbn_ (:, n_) = (OmgB \ Ckbn_) / fk(n_);
        Gan(:, :, n_) = I_N - fk(n_) / fz(n_) * (OmgA \ CZan_) / CKan_;
        Gbn(:, :, n_) = I_N - fz(n_) / fk(n_) * (OmgB \ CKbn_) / CZbn_;
        CKan(:, :, n_) = CKan_;
        CZbn(:, :, n_) = CZbn_;

    end

    for u_ = 1:n

        for n_ = 1:n
            Sun_ (:, u_, :, n_) = SFunction(u_, n_, Theta, OmgA, OmgB);
            SunGabn_(:, u_, :, n_) = SFunction(u_, n_, Theta, OmgA, OmgB) * [Gan(:, :, n_) zeros(N);zeros(N) Gbn(:, :, n_)];
        end

    end

    Yzan = reshape(Yzan_, [N * n 1]);
    Ykan = reshape(Ykan_, [N * n 1]);
    Ykbn = reshape(Ykbn_, [N * n 1]);
    Yzbn = reshape(Yzbn_, [N * n 1]);
    Sun = reshape(Sun_, [2 * N * n 2 * N * n]);
    SunGabn = reshape(SunGabn_, [2 * N * n 2 * N * n]);

    Yn = [Ykan; Yzbn] + Sun * ([Yzan; Ykbn] - [Ykan; Yzbn]);
    Xu_ = (I - SunGabn) \ Yn;

    XuN = reshape(Xu_, [N n 2]);
    XauN = XuN(:, :, 1)';
    XbuN = XuN(:, :, 2)';

    for n_ = 1:n
        a3nN(n_, :) = fk(n_) * (CKan(:, :, n_) \ XauN(n_, :)');
        b3nN(n_, :) = fz(n_) * (CZbn(:, :, n_) \ XbuN(n_, :)');
    end

end
