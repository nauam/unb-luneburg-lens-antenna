function index
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    diag =  "cart"; n = 100; N = 7; frq = 5*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    %Antena 1
    k0r = 20.94; tx = 0.990; ang1 = 4; ang2 = 2; r_s = tx; f_s = 1;
    [r__s,r,th_r,ep_r,mu_r] = centredProbeFed_discreteLuneburgLen(k0r,k0,r_s,tx,ang1,ang2,N);

    %Metodo de Regularização Analítica
    [~,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);

    %Diagrama de Campo Distante
    diag_f = farField(n,b3nN(:,1));
    diagram(diag,n,diag_f,"Diagrama de Campo Distante",35,30,5);
    
    %Diagrama de Diretividade
    diag_d = directivity(n,b3nN(:,1));
    diagram(diag,n,diag_d,"Diagrama de Diretividade",80,30,5);
    
end