function index
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    n = 120; N = 2; frq = 9*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    %Antena 0
    %k0r = 10; r1 = 1.000; r2 = 1.000; ang1 = 0; ang2 = 0; er = 1.00; r_s = r1; f_s = 1;

    %Antena 1
    %k0r =   k0; r1 = 0.038; r2 = 0.035; ang1 = 16; ang2 = 160; er = 1.23; r_s = r2; f_s = 1;

    %Antena 2
    k0r = 9.94; r1 = 1.000; r2 = 0.970; ang1 = 18; ang2 = 180; er = 1.30; r_s = r1; f_s = 1;

    [r__s,r,th_r,ep_r,mu_r] = centredProbeFed(k0r,k0,r_s,r1,r2,ang1,ang2,er);

    %Metodo de Regularização Analítica
    [~,~,a3nN,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
    
    %DIAGRAMA
    %diag =  "cart"; f_var1 = 12; f_var2 = 7; d_var1 = 12; d_var2 = 7;
    diag = "polar"; f_var1 = -40; f_var2 = 0; d_var1 = 0; d_var2 = 4;
    
    %Diagrama de Campo Distantes
    [diag_0,diag_pi_2] = farField(n,a3nN(:,1),b3nN(:,1));
    diagram(diag,n,diag_0,"Diagrama de Campo Distante",f_var1,f_var2);
    diagram(diag,n,diag_pi_2,"Diagrama de Campo Distante",f_var1,f_var2);
    
    %Diagrama de Diretividade
    diag_d = directivity(n,b3nN(:,1));
    diagram(diag,n,diag_d,"Diagrama de Diretividade",d_var1,d_var2);
    
end