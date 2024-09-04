function I=BTK_IE(E,VV, t, zplus, ND, mu_F, alpha1, alpha2,beta_BdG,kT,U)
q=1;
hbar=1;
beta_Dev = beta_BdG;
%% Surface Green's Functions
        ig1=(E+1j*zplus).*eye(2)-alpha1;        
        ig2=(E+1j*zplus).*eye(2)-alpha2;         
        
        gs1=inv(ig1);gs2=inv(ig2);

        change=1;
        %Recursive calculation of surface Green's function
        %tic
        while change >1e-8
        Gs=inv(ig1-beta_BdG'*gs1*beta_BdG);
        change=sum(sum(abs(Gs-gs1)))/(sum(sum(abs(gs1)+abs(Gs))));
        gs1=0.5*Gs+0.5*gs1;
        %gs1
        end
        %toc
        change=1;
        %tic
        while change >1e-8
        Gs=inv(ig2-beta_BdG*gs2*beta_BdG');
        change=sum(sum(abs(Gs-gs2)))/(sum(sum(abs(gs2)+abs(Gs))));
        gs2=0.5*Gs+0.5*gs2;
        %gs2
        end
        %toc                                              

%% Hamiltonian (for N-region)
    
    alpha = [2*t-mu_F+U 0; 0 -2*t+mu_F-U];
    H = zeros(2*ND, 2*ND);

    H(1,1) = alpha(1,1);
    H(1,2) = alpha(1,2);
    H(2,1) = alpha(2,1);
    H(2,2) = alpha(2,2);

    H(2*ND-1,2*ND-1) = alpha(1,1);
    H(2*ND-1,2*ND) = alpha(1,2);
    H(2*ND,2*ND-1) = alpha(2,1);
    H(2*ND,2*ND) = alpha(2,2);

    H(1,3) = beta_Dev(1,1);
    H(1,4) = beta_Dev(1,2);
    H(2,3) = beta_Dev(2,1);
    H(2,4) = beta_Dev(2,2);

    H(3,1) = beta_Dev(1,1);
    H(3,2) = beta_Dev(1,2);
    H(4,1) = beta_Dev(2,1);
    H(4,2) = beta_Dev(2,2);

    for n=3:2:2*ND-3                         
    	H(n,n)=alpha(1,1);
    	H(n,n+1)=alpha(1,2);
    	H(n+1,n)=alpha(2,1);
    	H(n+1,n+1)=alpha(2,2);

    	H(n,n+2) = beta_Dev(1,1);
    	H(n,n+3) = beta_Dev(1,2);
    	H(n+1,n+2) = beta_Dev(2,1);
    	H(n+1,n+3) = beta_Dev(2,2);

	H(n+2,n) = beta_Dev(1,1);
    	H(n+2,n+1) = beta_Dev(1,2);
    	H(n+3,n) = beta_Dev(2,1);
    	H(n+3,n+1) = beta_Dev(2,2);

    end 
    
 
%% Self-energy and broadening matrix 
   
    sig1=(beta_BdG'*gs1*beta_BdG);gam1=1i*(sig1-sig1');
    sig2=(beta_BdG*gs2*beta_BdG');gam2=1i*(sig2-sig2');

    Sigma1 = zeros(2*ND, 2*ND); Sigma2 = zeros(2*ND, 2*ND);
    Sigma1(1,1) = sig1(1,1); Sigma1(1,2) = sig1(1,2);
    Sigma1(2,1)=sig1(2,1); Sigma1(2,2)=sig1(2,2);
    Sigma2(2*ND-1,2*ND-1) = sig2(1,1); Sigma2(2*ND-1,2*ND) = sig2(1,2);
    Sigma2(2*ND,2*ND-1)=sig2(2,1); Sigma2(2*ND,2*ND)=sig2(2,2);

    Gamma1 = 1j*(Sigma1 - Sigma1');
    Gamma2 = 1j*(Sigma2 - Sigma2');
    
    GD =  inv((E+1j*zplus).*eye(2*ND) - H -1*real(Sigma1)- 1i*imag(Sigma1)-1*real(Sigma2) - 1i*imag(Sigma2));
    A = 1j*(GD-GD');
    DOS = trace(A);
    
    fermi1e = 1.0/(1.0 + exp((E-VV)/kT));
    fermi1h = 1.0/(1.0 + exp((E+VV)/kT));

    fermi2e = 1.0/(1.0 + exp(E/kT));
    fermi2h = 1.0/(1.0 + exp(E/kT));
    
    Fermi_matrix1 = [fermi1e 0;0 fermi1h];
    Fermi_matrix2 = [fermi2e 0;0 fermi2h];

    Fermi1 = kron(eye(ND),Fermi_matrix1);
    Fermi2 = kron(eye(ND),Fermi_matrix2);
    
    Sigma_corr = Gamma1*Fermi1+ Gamma2*Fermi2;
    G_corr = GD*Sigma_corr*GD';
    
    I_op = (1j*q/hbar)*(H*G_corr - G_corr*H);
    I = I_op(1,1)-I_op(2,2);
end