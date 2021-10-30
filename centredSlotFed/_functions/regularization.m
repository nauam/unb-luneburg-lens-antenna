function [Xau_,Xbu_,a3nN,b3nN] = regularization(n,N,k0,f_s,r__s,r,th_r,ep_r,mu_r)
double precision; format long;

    Ykan_   = zeros(N,n) ;  RunGan_ = zeros(N,n,N,n);   Gan  = zeros(N,N,n);    a3nN    = zeros(n,N);
    Yzan_   = zeros(N,n) ;  Run_    = zeros(N,n,N,n);   CKan = zeros(N,N,n);
    Ykbn_   = zeros(N,n) ;  SunGbn_ = zeros(N,n,N,n);   Gbn  = zeros(N,N,n);    b3nN    = zeros(n,N);
    Yzbn_   = zeros(N,n) ;  Sun_    = zeros(N,n,N,n);   CZbn = zeros(N,N,n);

    fk = ones(1,n);         I   = eye(N*n);     i0 = 1:N  ;     Ii = diag(ones(1,(N-1)),-1);
    fz = 1./(2.*(1:n)+1 );  I_N = eye(N)  ;     i1 = 2:N+1;     Is = diag(ones(1,(N-1)), 1);

    et_i0 = sqrt(mu_r(i0,i0)/ep_r(i0,i0));	k__i0 = k0*sqrt(mu_r(i0,i0)*ep_r(i0,i0));  d0 = double(f_s==i0-1)';
    et_i1 = sqrt(mu_r(i1,i1)/ep_r(i1,i1));  k__i1 = k0*sqrt(mu_r(i1,i1)*ep_r(i1,i1));  d1 = double(f_s==i1-1)'; 
    k___s = k__i1(f_s,f_s);                            

    [Z1_i0,Z3_i0,K1_i0,K3_i0] = sphericalFunction(k__i0*r   ,n);                Omg_a = - inv((inv(k__i0*et_i0) + inv(k__i1*et_i1))/(r*2)); 
    [Z1_i1,Z3_i1,K1_i1,K3_i1] = sphericalFunction(k__i1*r   ,n);                Omg_b = - inv((k__i0/et_i0 + k__i1/et_i1)*(r*2));
    [Z1__s,Z3__s,K1__s,K3__s] = sphericalFunction(k___s*r__s,n);

    for n_=1:n
        Z1Dan = Z1_i1(:,:,n_);          Z1Ean = Z1_i0(:,:,n_)*Ii;           Z1Fan = Z1_i0(:,:,n_);          Z1Tan = Z1Dan - Z1Ean;
        Z3Dan = Z3_i0(:,:,n_);          Z3Ean = Z3_i1(:,:,n_)*Is;           Z3Fan = Z3_i1(:,:,n_);          Z3Tan = Z3Dan - Z3Ean;
        K1Dan = K1_i1(:,:,n_)/et_i1;	K1Ean = K1_i0(:,:,n_)/et_i0*Ii;     K1Fan = K1_i0(:,:,n_)/et_i0;	K1Tan = K1Dan - K1Ean;
        K3Dan = K3_i0(:,:,n_)/et_i0;	K3Ean = K3_i1(:,:,n_)/et_i1*Is;     K3Fan = K3_i1(:,:,n_)/et_i1;	K3Tan = K3Dan - K3Ean;

        K1Dbn = K1_i1(:,:,n_);          K1Ebn = K1_i0(:,:,n_)*Ii;           K1Fbn = K1_i0(:,:,n_);          K1Tbn = K1Dbn - K1Ebn;
        K3Dbn = K3_i0(:,:,n_);          K3Ebn = K3_i1(:,:,n_)*Is;           K3Fbn = K3_i1(:,:,n_);          K3Tbn = K3Dbn - K3Ebn;
        Z1Dbn = Z1_i1(:,:,n_)/et_i1;	Z1Ebn = Z1_i0(:,:,n_)/et_i0*Ii;     Z1Fbn = Z1_i0(:,:,n_)/et_i0;	Z1Tbn = Z1Dbn - Z1Ebn;
        Z3Dbn = Z3_i0(:,:,n_)/et_i0;	Z3Ebn = Z3_i1(:,:,n_)/et_i1*Is;     Z3Fbn = Z3_i1(:,:,n_)/et_i1;	Z3Tbn = Z3Dbn - Z3Ebn;

        CKan_ = K3Tan - K1Tan/Z1Tan*Z3Tan;   A1s =  K1__s(1,1,n_)*d1;       CKan(:,:,n_) = CKan_;
        CZan_ = Z3Dan + Z1Ean/Z1Tan*Z3Tan;   A3s =  K3__s(1,1,n_)*d0;       Gan(:,:,n_) = I_N - fk(n_)/fz(n_)*(Omg_a\CZan_)/CKan_;

        CZbn_ = Z3Tbn - Z1Tbn/K1Tbn*K3Tbn;   B1s =  Z1__s(1,1,n_)*d1;       CZbn(:,:,n_) = CZbn_;
        CKbn_ = K3Dbn + K1Ebn/K1Tbn*K3Tbn;   B3s =  Z3__s(1,1,n_)*d0;       Gbn(:,:,n_) = I_N - fz(n_)/fk(n_)*(Omg_b\CKbn_)/CZbn_;

        Czan_ =          Z1Ean/Z1Tan*Z3Fan *A1s - (eye(N) + Z1Ean/Z1Tan)*Z1Fan *A3s;  Yzan_ (:,n_) = (Omg_a\Czan_)/fz(n_);
        Ckan_ = (K3Fan - K1Tan/Z1Tan*Z3Fan)*A1s - ( K1Fan - K1Tan/Z1Tan *Z1Fan)*A3s;  Ykan_ (:,n_) =       Ckan_/fk(n_);
    
        Ckbn_ =          K1Ebn/K1Tbn*K3Fbn *B1s - (eye(N) + K1Ebn/K1Tbn)*K1Fbn *B3s;  Ykbn_ (:,n_) = (Omg_b\Ckbn_)/fk(n_);
        Czbn_ = (Z3Fbn - Z1Tbn/K1Tbn*K3Fbn)*B1s - ( Z1Fbn - Z1Tbn/K1Tbn *K1Fbn)*B3s;  Yzbn_ (:,n_) =       Czbn_/fz(n_);
    end

    for u_=1:n
        for n_=1:n
           Run_   (:,u_,:,n_) = func_r(u_,n_,th_r);
           RunGan_(:,u_,:,n_) = func_r(u_,n_,th_r)*Gan(:,:,n_);
           Sun_   (:,u_,:,n_) = func_s(u_,n_,th_r);
           SunGbn_(:,u_,:,n_) = func_s(u_,n_,th_r)*Gbn(:,:,n_);
        end
    end
    
    Run    = reshape(Run_   ,[N*n N*n]);        Yzan = reshape(Yzan_,[N*n 1]);          Ykan = reshape(Ykan_,[N*n 1]);
    RunGan = reshape(RunGan_,[N*n N*n]);        Yan  = Ykan + Run*(Yzan - Ykan);        
    Xau_    = (I-RunGan)\Yan;                    Xau   = reshape(Xau_,[N n]);

    Sun    = reshape(Sun_   ,[N*n N*n]);        Ykbn = reshape(Ykbn_,[N*n 1]);          Yzbn = reshape(Yzbn_,[N*n 1]); 
    SunGbn = reshape(SunGbn_,[N*n N*n]);        Ybn  = Yzbn + Sun*(Ykbn - Yzbn);        
    Xbu_    = (I-SunGbn)\Ybn;                    Xbu   = reshape(Xbu_,[N n]);

    for n_ = 1:n
        a3nN(n_,:) = fk(n_)*(CKan(:,:,n_)\Xau(:,n_));
        b3nN(n_,:) = fz(n_)*(CZbn(:,:,n_)\Xbu(:,n_));
    end
    warning('query','all');
end