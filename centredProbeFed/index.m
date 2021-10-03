function index
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    diag =  "cart"; n = 20; N = 2; frq = 9*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    %Antena 1
    %k0r =   k0; r1 = 0.038; r2 = 0.035; ang1 = 16; ang2 = 160; er = 1.23; r_s = r2; f_s = 1;

    %Antena 2
    k0r = 9.95; r1 = 1.000; r2 = 0.970; ang1 = 18; ang2 = 180; er = 1.30; r_s = r2; f_s = 1;

    [r__s,r,th_r,ep_r,mu_r] = centredProbeFed(k0r,k0,r_s,r1,r2,ang1,ang2,er);

    %Metodo de Regularização Analítica
    [~,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);

    %Diagrama de Campo Distante
    diag_f = farField(n,b3nN(:,1));
    diagram(diag,n,diag_f);
    
    %Diagrama de Diretividade
    diag_d = directivity(n,b3nN(:,1));
    diagram(diag,n,diag_d);
    
end