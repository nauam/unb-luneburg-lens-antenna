function an_directivity
clear all; double precision; format long; addpath(genpath('_functions'))
    n = 50; N = 2; frq = 9*10^9; k0  = 2*pi*frq/299792458;
    
    %Estrutura da Antena
    r1 = 1.000; r2 = 0.970; ang1 = 18; ang2 = 180; er = 1.30; r_s = r2; f_s = 1;
    k0r_max = 70;
    k0r = 0:0.01:k0r_max;
    Pr_norm = zeros(length(k0r),1);                                         i_ = 1;
    
    f = waitbar(0,strcat('k0r = ',num2str(length(k0r))));
    for k0r_ = k0r
        [r__s,r,th_r,ep_r,mu_r] = str_probeFed(k0r_,k0,r_s,r1,r2,ang1,ang2,er);

        %Metodo de Regularização Analítica
        [~,b3nN] = def_regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
        
        %Potência de radiação
        Pr_norm(i_,:) = def_radiatedPower(n,b3nN(:,1));                     i_ = i_ + 1;
        
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
    
    %Gráfico
    figure(1)
    semilogy(1:length(k0r),Pr_norm,'Linewidth',2);  hold on
    findpeaks(-real(Pr_norm),'MinPeakProminence',1)
    findpeaks( real(Pr_norm),'MinPeakProminence',1)
    
end

