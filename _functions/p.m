function P = p(u_,n_,theta)
double precision; format long; 

    P = theta/2*(sinc((u_-n_)*theta/pi) - sinc((u_+n_+1)*theta/pi));
    
end