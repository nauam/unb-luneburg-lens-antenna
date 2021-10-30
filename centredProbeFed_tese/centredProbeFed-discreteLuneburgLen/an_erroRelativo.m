function an_erroRelativo
    clear all; double precision; format long; 
    
    %Estrutura da Antena
    n = 80; N = 7; k0r = [20.94]; frq = 5*10^9; k0  = 2*pi*frq/299792458; X = zeros((n+1)*N,n+1); erro = zeros(n,length(k0r)); 
    
                ang1 = 4; f_s = 1; 
    tx = 0.97; ang2 = 2; r_s = tx; 
    
    i_ = 1;
    
    for k0r_ = k0r
        
        [r__s,r,th_r,ep_r,mu_r] = centredProbeFed_discreteLuneburgLen(k0r_,k0,r_s,tx,ang1,ang2,N);

        %Metodo de Regulariza��o Anal�tica
        f = waitbar(0,strcat('n = ',num2str(n)));
        for n_ = 1:n+1
            [X(1:n_*N,n_),~] = regularization(n_,N,k0,f_s,r__s,r,th_r,ep_r,mu_r);
            waitbar(n_/n, f, sprintf(strcat('n_ = �',num2str(n_),'` �',num2str(n+1),'` e k0r = �',num2str(k0r_),'�')))
        end
        delete(f)

        %Erro Relativo
        for n_ = 1:n
            erro(n_,i_) = norm(X(1:n_*N,n_+1)-X(1:n_*N,n_),inf)./norm(X(1:n_*N,n_+1),inf);
        end

        legend_str{i_} = strcat('k0r = ',num2str(k0r_)); 
        i_ = i_+1;
    end
    
    %Gr�fico
    figure(3)
    semilogy(1:n,erro,'Linewidth',2);  hold on
    grid on; legend(legend_str);
    title(strcat('Lente de Luneburg Discreta,  = 18�,  r_2/r_1 = 0.97 e Th_2   = ',num2str(ang2))); 
    xlabel('Tuncamento n'); 
    ylabel('Erro Relativo e_r(n)');
end