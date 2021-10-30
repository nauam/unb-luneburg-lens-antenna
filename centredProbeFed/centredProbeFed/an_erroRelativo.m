function an_erroRelativo
    clear all; double precision; format long; 
    addpath(genpath('_functions')); addpath(genpath('_structure')); addpath(genpath('_diagram'));
    
    %Estrutura da Antena
    n = 80; N = 2; k0r = [9.95]; frq = 9*10^9; k0  = 2*pi*frq/299792458; X = zeros((n+1)*N,n+1); erro = zeros(n,length(k0r)); 
    
    r1 = 1.000; ang1 = 18;  f_s = 1; er = 1.30; 
    r2 = 0.970; ang2 = 180; r_s = r2; 
    
    i_ = 1;
    
    for k0r_ = k0r
        
        [r__s,r,th_r,ep_r,mu_r] = centredProbeFed(k0r_,k0,r_s,r1,r2,ang1,ang2,er);

        %Metodo de Regularização Analítica
        f = waitbar(0,strcat('n = ',num2str(n)));
        for n_ = 1:n+1
            [X(1:n_*N,n_),~] = regularization(n_,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
            waitbar(n_/n, f, sprintf(strcat('n_ = ´',num2str(n_),'` ´',num2str(n+1),'` e k0r = ´',num2str(k0r_),'´')))
        end
        delete(f)

        %Erro Relativo
        for n_ = 1:n
            erro(n_,i_) = norm(X(1:n_*N,n_+1)-X(1:n_*N,n_),1)./norm(X(1:n_*N,n_+1),1);
        end

        legend_str{i_} = strcat('k0r = ',num2str(k0r_)); 
        i_ = i_+1;
    end
    
    %Gráfico
    figure(1)
    semilogy(1:n,erro,'Linewidth',2);  hold on
    grid on; legend(legend_str);
    title(strcat('Lente de Luneburg Discreta,  = 18°,  r_2/r_1 = 0.97 e Th_2   = ',num2str(ang2))); 
    xlabel('Tuncamento n'); 
    ylabel('Erro Relativo e_r(n)');
end
