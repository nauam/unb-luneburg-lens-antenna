function Index
    addpath(genpath('_function')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    clear variables; double precision; format long;
    
    n = 120;
    
    %----------------------------(SELECT STRUCTURE AND RUN)----------------------------------------
    %----------------------------------Antenna Structure-------------------------------------------
        %%------------------------------Dipolo magétrico--------------------------------------------
        %%%SETTINGS <=> [NLENS] = [0]
        %antenna = "Dipolo"; settings = 1; %%%choose settings = [1]

        %%--------------------------------LLD + dipolo---------------------------------------------
        %%%SETTINGS <=> [NLENS] = [1, 2, 3, 20]
        antenna = "Dipolo"; settings = 1; %%%choose settings = [2; 3; 4; 5]
        
        %%--------------------------Antenna + LLD + dipolo-----------------------------------------
        %%%SETTINGS <=> [NLENS, K0R] = [[0, 11.05]; [1, 10.65]; [2, 10.70]; [3, 10.85]; [20, 11.00]]
        %antenna = "Antenna_LLD_Dipolo"; settings = 1; %%%choose settings = [1; 2; 3; 4; 5]

        %%---------------------Analyze - Antenna + LLD + dipolo------------------------------------
        %antenna = "Analyze_Antenna_LLD_Dipolo"; settings = [1; 2; 3; 4; 5];
    %----------------------------------------------------------------------------------------------

    if antenna == "Analyze_Antenna_LLD_Dipolo"
        %Relative Error
        RelativeErrorDiagram(n, 7e6 + n, antenna, [1, 2, 3, 4, 5], 'N');
        
        %Radiated Power
        RadiatedPowerDiagram(n, 8e6 + n, antenna, [1, 2, 3, 4, 5], 18, 0.05);
    else
        [Nlens, k0r, omega, r, ang, epsilonR, source, AxisF] = Structure(antenna, settings);
        [N, R, Theta, Epsilon, Mu] = LuneburgLensStructure(Nlens, k0r, omega, r, ang, epsilonR);

        %Method Analytical Regularization
        [~, ~, a3nN, b3nN, Gan, Gbn] = RegularizationFunction(n, N, omega, source, R, Theta, Epsilon, Mu);

        %Singularity
        titleFigure = "Singularidade A";
        Singularity(1e5 + n, titleFigure, Gan, n);

        %Singularity
        titleFigure = "Singularidade B";
        Singularity(2e5 + n, titleFigure, Gbn, n);

        %Diagram
        m = 1; a3e0mn = a3nN(:, 1); a3o0mn = 0; b3e0mn = b3nN(:, 1); b3o0mn = 0;

        %Far Field Diagram
        [dTh, dPh] = FarFieldFunction(m, n, a3e0mn, a3o0mn, b3e0mn, b3o0mn, pi/2);

        Diagram(3e6 + n, dTh, "Diagrama de Campo Distante - Componente Principal - HTheta", AxisF, 'b');
        Diagram(4e6 + n, dPh, "Diagrama de Campo Distante - Componente Principal - HTheta", AxisF, 'g');

        [dTh, dPh] = FarFieldFunction(m, n, a3e0mn, a3o0mn, b3e0mn, b3o0mn, 0);

        Diagram(5e6 + n, dTh, "Diagrama de Campo Distante - Componente Cruzada - HTheta", AxisF, 'b');
        Diagram(6e6 + n, dPh, "Diagrama de Campo Distante - Componente Cruzada - HTheta", AxisF, 'g');
    end
end
