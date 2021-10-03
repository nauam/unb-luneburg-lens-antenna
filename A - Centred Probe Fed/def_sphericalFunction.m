function [Z1,Z3,K1,K3] = def_sphericalFunction(z,n)
double precision; format long;

N = length(z);  Z1 = zeros(N,N,n);	Z3 = zeros(N,N,n);	K1 = zeros(N,N,n);	K3 = zeros(N,N,n);

n_     = (1:n); 
n_or05 = (1:(n+1)) + 1/2;
    
for i_ = 1:N
    z_i = z(i_,i_);
    
    jn = besselj(n_or05,z_i).*sqrt(pi/(2*z_i));
    yn = bessely(n_or05,z_i).*sqrt(pi/(2*z_i));
    hn = jn + 1i*yn;

    Z1(i_,i_,:) = jn(n_);
    Z3(i_,i_,:) = hn(n_);
    K1(i_,i_,:) = jn(n_).*(n_+1)/z_i - jn(n_+1);
    K3(i_,i_,:) = hn(n_).*(n_+1)/z_i - hn(n_+1);
end
end

