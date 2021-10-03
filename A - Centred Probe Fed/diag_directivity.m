function diag_D = diag_directivity(n,b3n1)
double precision; format long;

    %Campo distante
    Ff = def_farField(n,b3n1);
    
    %Potência de radiação
    Pr_norm = def_radiatedPower(n,b3n1);
    
    %Diagrama de diretividade
    diag_D = Ff.*conj(Ff)/Pr_norm;

end