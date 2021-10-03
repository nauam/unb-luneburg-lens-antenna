function Q = def_q(u_,n_,theta)
double precision; format long; 

    Q = theta/2*(sinc((u_-n_)*theta/pi) + sinc((u_+n_+1)*theta/pi));
    
end