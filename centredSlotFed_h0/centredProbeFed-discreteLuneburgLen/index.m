function index
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    diag =  "polar"; n = 100; frq = 5*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    %Antena 1
    %N = 2; k0r = 10; tx = 1.000; ang1 = 0; ang2 = 0; r_s = 1; f_s = 1;
    
    %Antena 2
    N = 20; k0r = 20.94; tx = 0.990; ang1 = 4; ang2 = 2; r_s = 1; f_s = 1;
    
    [r__s,r,th_r,ep_r,mu_r] = centredProbeFed_discreteLuneburgLen(k0r,k0,r_s,tx,ang1,ang2,N);

    %Metodo de Regularização Analítica
    [~,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
    
    %DIAGRAMA
    %diag =  "cart"; f_var1 = 12; f_var2 = 11; d_var1 = 12; d_var2 = 8;
    diag = "polar"; f_var1 = -60; f_var2 = 0; d_var1 = 0; d_var2 = 40;

    %Diagrama de Campo Distante
    diag_f = farField(n,b3nN(:,1));
    diagram(diag,n,diag_f,"Diagrama de Campo Distante",f_var1,f_var2);
    
    %Diagrama de Diretividade
    diag_d = directivity(n,b3nN(:,1));
    diagram(diag,n,diag_d,"Diagrama de Diretividade",d_var1,d_var2);
    
end