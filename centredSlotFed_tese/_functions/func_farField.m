function [Ff_0,Ff_pi_2] = func_farField(n,a3n1,b3n1)
double precision; format long;

    P1n = zeros(1001,n);
    dP1n = zeros(1001,n);
    
    cos_t = cos(0:pi/500:2*pi);
    sin_t = sin(0:pi/500:2*pi);
    
    for n_ = 1:n
         Pmn = legendre           (n_,cos_t);  P1n(:,n_) =    Pmn(2,:)./sin_t;  % Pmn(cos(t))/sin(t)
        dPmn = legendre_derivative(n_,cos_t); dP1n(:,n_) = - dPmn(2,:).*sin_t;  %dPmn(cos(t))/dt
    end
    
    i2np_nnp = diag((-1i).^((1:n)+1));
    
    Ff_0     = 1i*dP1n*i2np_nnp*a3n1 +    P1n*i2np_nnp*b3n1;
    Ff_pi_2  =    dP1n*i2np_nnp*b3n1 + 1i*P1n*i2np_nnp*a3n1;
    
end