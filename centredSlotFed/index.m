function index
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    dgr =  "cart"; n = 120; N = 2; frq = 9*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    %Antena 1
    %k0r =   k0; r1 = 0.038; r2 = 0.035; ang1 = 16; ang2 = 160; er = 1.23; r_s = r2; f_s = 1;

    %Antena 2
    k0r = 5.13; r1 = 1.000; r2 = 0.990; ang1 = 18; ang2 = 180; er = 1.30; r_s = r2; f_s = 1;

    [r__s,r,th_r,ep_r,mu_r] = centredSlotFed(k0r,k0,r_s,r1,r2,ang1,ang2,er);

    %Metodo de Regularização Analítica
    [~,~,a3nN,b3nN] = regularizationSlot(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);

    %Diagrama de Campo Distante
    diag_f = farField(n,b3nN(:,1));
    diagram(dgr,n,diag_f,"Diagrama de Campo Distante",35,30,5);
    
    %Diagrama de Diretividade
    diag_d = directivity(n,b3nN(:,1));
    diagram(dgr,n,diag_d,"Diagrama de Diretividade",80,30,5);
    
end