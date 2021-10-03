function diag_Ff = diag_farField(n,b3n1)
double precision; format long;
    
    %Campo distante
    Ff = def_farField(n,b3n1);
    
    %Diagrama de campo distante
    diag_Ff = abs(Ff)/max(abs(Ff));
    
end