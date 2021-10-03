function Ff = def_farField(n,b3n1)
double precision; format long;

    P1n = zeros(1001,n);
    for n_ = 1:n
        Pmn = legendre(n_,cos(0:pi/500:2*pi));  P1n(:,n_) = Pmn(2,:);
    end
    Ff  = (P1n*diag((-1i).^(1:n).*(2*(1:n)+1))*b3n1);
    
end