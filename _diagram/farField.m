function Ff = farField(n,b3n1)
double precision; format long; addpath(genpath('_functions'))
    
    %Campo distante
    func_Ff = func_farField(n,b3n1);
    
    %Diagrama de campo distante
    Ff = abs(func_Ff)/max(abs(func_Ff));
    
end