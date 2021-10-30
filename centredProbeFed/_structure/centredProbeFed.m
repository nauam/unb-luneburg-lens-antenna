function [r__s,r,th_r,ep_r,mu_r] = centredProbeFed(k0r,k0,r_s,r1,r2,ang1,ang2,er)
double precision; format long;

    r__s = r_s*k0r/k0;
    r    = diag([             r1       r2])*k0r/k0;
    th_r = diag([           ang1     ang2])*pi/180;
    ep_r = diag([  1.000      er    1.000]);
    mu_r = diag([  1.000   1.000    1.000]);

end