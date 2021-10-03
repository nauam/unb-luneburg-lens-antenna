function Pr_norm = def_radiatedPower(n,b3n1)
double precision; format long;

    Pr_norm = ((1:n).*((1:n)+1).*(2*(1:n)+1))*(b3n1.*conj(b3n1));
    
end