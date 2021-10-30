function [r__s,r,th_r,ep_r,mu_r] = centredProbeFed_discreteLuneburgLen(k0r,k0,r_s,tx,ang1,ang2,N)
    double precision; format long;

       i = (N-1:-1:1)';
    r__i = i/(N-1)*tx;
    th_i = zeros((N-2),1);
    ep_i = 2 - ((2*i-1)/(2*(N-1))).^2;
    mu_i = ones(N-1,1);

    r__s = r_s*k0r/k0;
    r    = diag([       1.000; r__i])*k0r/k0;
    th_r = diag([ ang1;  ang2; th_i])*pi/180;   
    ep_r = diag([1.000; 1.000; ep_i])       ;
    mu_r = diag([1.000; 1.000; mu_i])       ;

end