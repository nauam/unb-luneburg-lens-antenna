function an_radiatedPowerAndDirectivity
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    
    %Estrutura da Antena
    frq      = 9*10^9;          k0 = 2*pi*frq/299792458;        n   = 120;  N  =    2;
    k0r_max  =    100;          r1 = 1.000;     ang1 =  18;     f_s =  1;   er = 1.30; 
    k0r_step =   0.01;          r2 = 0.970;     ang2 = 180;     r_s = r2; 
    k0r = 0:k0r_step:k0r_max;  
    
    Pr_norm = zeros(length(k0r),1);
    direct  = zeros(length(k0r),1);             i_ = 1;                
    
    f = waitbar(0,strcat('k0r = ',num2str(length(k0r))));
    for k0r_ = k0r
        [r__s,r,th_r,ep_r,mu_r] = centredProbeFed(k0r_,k0,r_s,r1,r2,ang1,ang2,er);

        %Metodo de Regularização Analítica
        [~,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
        
        %Potência de radiação
        Pr_norm(i_,:) = func_radiatedPower(n,b3nN(:,1));
        
        %Maxima diretividade
        direct(i_,:) = max(directivity(n,b3nN(:,1)));                     i_ = i_ + 1;
        
        %Barra de carregamento
        waitbar(k0r_/k0r_max, f, sprintf(strcat('k0r = ´',num2str(k0r_),'` ´',num2str(k0r_max),'´')))
    end
    delete(f)
    
    %Picos e Vales
    [v_value,v_key] = findpeaks(-real(Pr_norm),'MinPeakProminence',1);
    [p_value,p_key] = findpeaks( real(Pr_norm),'MinPeakProminence',1);
    
    disp(v_value) 
    disp(v_key/100)
    disp(p_value)
    disp(p_key/100)
    
    %Gráfico Potência de radiação
    figure(1)
    semilogy(1:length(k0r),Pr_norm,'Linewidth',2);  hold on
    findpeaks(-real(Pr_norm),'MinPeakProminence',1)
    findpeaks( real(Pr_norm),'MinPeakProminence',1)
    
    %Gráfico Maxima diretividade
    figure(2)
    semilogy(1:length(k0r),direct,'Linewidth',2);  hold on
    findpeaks(-real(direct),'MinPeakProminence',1)
    findpeaks( real(direct),'MinPeakProminence',1)
    
end

