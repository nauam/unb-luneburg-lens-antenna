function [Nlens, k0r, omega, r, ang, epsilonR, source, AxisF] = Structure(antenna, settings)
    double precision;
    format long;

    NLENS = [0, 1, 2, 3, 20];
    K0R = [11.05, 10.65, 10.70, 10.85,11.00];
    
    if antenna == "Dipolo" %Lente de Lüneburg discreta associada a um dipolo elétrico
        Nlens = NLENS(settings);
        k0r = 10.00;
        ang = [0, 0];
        tx = 1;
    else%if antenna == "Antenna_LLD_Dipolo" %theseSeb p159
        Nlens = NLENS(settings);
        k0r = K0R(settings);
        ang = [36, 18];
        tx = 0.97;
    end
    
    frq = 5 * 10^9;
    omega = 2 * pi * frq;

    rs = 2;
    r = [1 / tx, 1];
    source = [1, rs];
    epsilonR = 1.00;
    
    maxF = 0;
    stepF = 5;
    minF = -40;
    diagram = "polar";
    normalization = "dB";
    AxisF = [normalization, diagram, minF, maxF, stepF];

end
