function [N, R, Theta, Epsilon, Mu] = LuneburgLensStructure(Nlens, k0r, omega, r, ang, epsilonR)
    double precision;
    format long;

    c = 299792458;
    k0 = omega / c;
    mu0 = 4 * pi * 1e-7;
    epsilon0 = 1 / (mu0 * c^2);

    N = 2;
    lensR = r(2);
    lensEpsilon = 1.000;

    if Nlens > 0
        N = 1 + Nlens;
        lensI = (Nlens:-1:1)';
        lensR = lensI / Nlens;
        lensEpsilon = 2 - ((2 * lensI - 1) / (2 * Nlens)).^2;
    end

    R = diag([r(1); lensR]) * k0r / k0;
    Theta = diag([ang(1); ang(2); zeros((Nlens - 1), 1)]) * pi / 180;
    Epsilon = diag([1.000; epsilonR; lensEpsilon]) * epsilon0;
    Mu = diag([1.000; 1.000; 1.000; ones(Nlens - 1, 1)]) * mu0;

end
