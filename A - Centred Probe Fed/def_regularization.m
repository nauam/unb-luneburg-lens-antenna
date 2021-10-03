function [Xu_,b3nN] = def_regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r)
double precision; format long;

    Ykbn_   = zeros(N,n) ;  SunGbn_ = zeros(N,n,N,n);   Gbn  = zeros(N,N,n);    b3nN    = zeros(n,N);
    Yzbn_   = zeros(N,n) ;  Sun_    = zeros(N,n,N,n);   CZbn = zeros(N,N,n);

    fk = 1./(      2*(1:n)+1 );             I   = eye(N*n);     i0 = 1:N  ;     Ii = diag(ones(1,(N-1)),-1);
    fz = 1./((1:n).*((1:n)+1));             I_N = eye(N)  ;     i1 = 2:N+1;     Is = diag(ones(1,(N-1)), 1);

    et_i0 = sqrt(mu_r(i0,i0)/ep_r(i0,i0));	k__i0 = k0*sqrt(mu_r(i0,i0)*ep_r(i0,i0));  d0 = double(f_s==i0-1)';
    et_i1 = sqrt(mu_r(i1,i1)/ep_r(i1,i1));  k__i1 = k0*sqrt(mu_r(i1,i1)*ep_r(i1,i1));  d1 = double(f_s==i1-1)'; k___s = k__i1(f_s,f_s);                            

    [Z1_i0,Z3_i0,K1_i0,K3_i0] = def_sphericalFunction(k__i0*r   ,n);                 Omg_b = - inv(r/2*(k__i0/et_i0 + k__i1/et_i1));
    [Z1_i1,Z3_i1,K1_i1,K3_i1] = def_sphericalFunction(k__i1*r   ,n);
    [Z1__s,Z3__s,~    ,~    ] = def_sphericalFunction(k___s*r__s,n);

    for n_=1:n
        K1Dbn = K1_i1(:,:,n_);          K1Ebn = K1_i0(:,:,n_)*Ii;           K1Fbn = K1_i0(:,:,n_);          K1Tbn = K1Dbn - K1Ebn;
        K3Dbn = K3_i0(:,:,n_);          K3Ebn = K3_i1(:,:,n_)*Is;           K3Fbn = K3_i1(:,:,n_);          K3Tbn = K3Dbn - K3Ebn;
        Z1Dbn = Z1_i1(:,:,n_)/et_i1;	Z1Ebn = Z1_i0(:,:,n_)/et_i0*Ii;     Z1Fbn = Z1_i0(:,:,n_)/et_i1;	Z1Tbn = Z1Dbn - Z1Ebn;
        Z3Dbn = Z3_i0(:,:,n_)/et_i0;	Z3Ebn = Z3_i1(:,:,n_)/et_i1*Is;     Z3Fbn = Z3_i1(:,:,n_)/et_i0;	Z3Tbn = Z3Dbn - Z3Ebn;

        CZbn_ = Z3Tbn - Z1Tbn/K1Tbn*K3Tbn;   B1s =  Z1__s(1,1,n_)*d1;       CZbn(:,:,n_) = CZbn_;
        CKbn_ = K3Dbn + K1Ebn/K1Tbn*K3Tbn;   B3s =  Z3__s(1,1,n_)*d0;       Gbn(:,:,n_) = I_N - fz(n_)/fk(n_)*(Omg_b\CKbn_)/CZbn_;

        Ckbn_ =          K1Ebn/K1Tbn*K3Fbn *B1s - (eye(N) + K1Ebn/K1Tbn)*K1Fbn *B3s;  Ykbn_ (:,n_) = (Omg_b\Ckbn_)/fk(n_);
        Czbn_ = (Z3Fbn - Z1Tbn/K1Tbn*K3Fbn)*B1s - ( Z1Fbn - Z1Tbn/K1Tbn *K1Fbn)*B3s;  Yzbn_ (:,n_) =       Czbn_/fz(n_);
    end

    for u_=1:n
        for n_=1:n
           Sun_   (:,u_,:,n_) = def_s(u_,n_,th_r);
           SunGbn_(:,u_,:,n_) = def_s(u_,n_,th_r)*Gbn(:,:,n_);
        end
    end

    Sun    = reshape(Sun_   ,[N*n N*n]);        Ykbn = reshape(Ykbn_,[N*n 1]);          Yzbn = reshape(Yzbn_,[N*n 1]); 
    SunGbn = reshape(SunGbn_,[N*n N*n]);        Ybn  = Yzbn + Sun*(Ykbn - Yzbn);        
    Xu_    = (I-SunGbn)\Ybn;                    Xu   = reshape(Xu_,[N n]);

    for n_ = 1:n
        b3nN(n_,:) = fz(n_)*(CZbn(:,:,n_)\Xu(:,n_));
    end
    warning('query','all');
end