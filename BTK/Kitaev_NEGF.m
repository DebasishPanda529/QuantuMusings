clear all
hbar=1;
q=1;
kT = 0.0000001;
Delta1=0;
Delta2=1e-3;
t=1.7;
mu_F=2*t;
UB=0; % zero barrier

E=linspace(-2*Delta2, 2*Delta2, 1001);
zplus=1e-6;
ND=2;
V = linspace(0,0.2*Delta2,11);
j=1;


for k=1:length(V)% = [0:dV:2*Delta2]          

  for i=1:length(E)
      %VV = 1*Delta2/2;
  tic

    alpha1 = [2*t-mu_F(j) Delta1; Delta1' -2*t+mu_F(j)];
    alpha2=[2*t-mu_F(j) Delta2;Delta2' -2*t+mu_F(j)];
    beta_BdG=[-t 0;0 t];


            U=0;
 

             I_E =  @(E) BTK_IE(E,V(k), t, zplus, ND, mu_F(j), alpha1, alpha2,beta_BdG,kT,U);
             %I_E(i) =  BTK_IE(E(i),V, t, zplus, ND, mu_F, alpha1, alpha2,beta_BdG,kT,U);
             I1(k) = integral(I_E, -1*abs(V(k)), 1*abs(V(k)),'AbsTol',1e-8,'ArrayValued',true);
 
  i

  toc

  end
end
figure(1)
plot(V./Delta2,real(I1)./Delta2,'r','linewidth',4)
title('I-V Characteristics','fontSize',40,'interpreter','latex')
ylabel('Current$(e\Delta_0/\hbar)$','interpreter','latex','fontsize',40)
xlabel('$eV/\Delta_0$','fontSize',40,'interpreter','latex')
%legend('U=0', 'U=0.5t', 'U=t', 'U=2t')
set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
  %%

 figure(2)
 hold on
dV = V(2)-V(1);
G1=diff(real(I1))./dV;

plot(V(1:length(V)-1)./Delta2,G1,'b','linewidth',4)
title('Conductance','fontSize',40,'interpreter','latex')
ylabel('G$(e^2/\hbar)$','interpreter','latex','fontsize',40)
xlabel('$eV/\Delta_0$','fontSize',40,'interpreter','latex')
%legend('U=0', 'U=0.5t', 'U=t', 'U=2t')
set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
