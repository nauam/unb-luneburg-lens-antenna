function [Ff_0,Ff_pi_2] = farField(n,a3n1,b3n1)
double precision; format long; addpath(genpath('_functions'))
    
    %Campo distante
    [func_Ff_0,func_Ff_pi_2] = func_farField(n,a3n1,b3n1);
    
    %Diagrama de campo distante
    Ff_0    = 20*log10(abs(func_Ff_0   )/max(abs(func_Ff_0   )));
    Ff_pi_2 = 20*log10(abs(func_Ff_pi_2)/max(abs(func_Ff_pi_2)));
    
    Ff_0(isinf(Ff_0)) = NaN;
    Ff_pi_2(isinf(Ff_pi_2)) = NaN;
    
end