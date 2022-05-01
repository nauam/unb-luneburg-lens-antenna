function S = SFunction(u_, n_, Theta, OmgA, OmgB)
    double precision;
    format long;

    s11 = zeros(length(Theta));
    s12 = zeros(length(Theta));
    s21 = zeros(length(Theta));
    s22 = zeros(length(Theta));

    for i = 1:length(Theta)
        theta = Theta(i, i);
        omgA = OmgA(i, i);
        omgB = OmgB(i, i);
        
        if theta == 0
            s = [0 0;0 0];
        elseif theta == pi
            if u_ ==  n_
                s = [1 0;0 1];
            else
                s = [0 0;0 0];
            end
        else
            sun = [PFunction(u_, n_, theta) 0 ; 0 QFunction(u_, n_, theta)];
            su0 = [PFunction(u_, 0, theta)/omgA -PFunction(u_, 0, theta); QFunction(u_, 0, theta)/omgB QFunction(u_, 0, theta)];
            s00 = [PFunction(0, 0, theta)/omgA PFunction(0, 0, pi)-PFunction(0, 0, theta); QFunction(0, 0, theta)/omgB QFunction(0, 0, theta)-QFunction(0, 0, pi)];
            s0n = [PFunction(0, n_, theta) 0 ; 0 QFunction(0, n_, theta)];
            
            s = 2 / pi * (sun - (su0 / s00) * s0n);
        end
        
        s11(i,i) = s(1, 1);
        s12(i,i) = s(1, 2);
        s21(i,i) = s(2, 1);
        s22(i,i) = s(2, 2);

    end
    
    S = [s11 s12; s21 s22];
    
end
