clear all
hbar=1;
q=1;
kT = 0.0000001;
Delta1=0;
Delta2=1e-3;
t=1.7;
mu_F=2*t;

E=linspace(-2*Delta2, 2*Delta2, 1001);
zplus=1e-6;
ND=2;
V = linspace(0,2*Delta2,61);    
j=1;

UB_array = [0, 0.5*t, t, 2*t];
I1_all = zeros(length(UB_array), length(V));
G1_all = zeros(length(UB_array), length(V)-1);

colors = ['r', 'g', 'b', 'k'];
labels = {'1', '2', '3', '4'};

for u=1:length(UB_array)
    UB = UB_array(u);

    for k=1:length(V)% = [0:dV:2*Delta2]          

      for i=1:length(E)
          %VV = 1*Delta2/2;
      tic
    
        alpha1 = [2*t-mu_F(j) Delta1; Delta1' -2*t+mu_F(j)];
        alpha2 = [2*t-mu_F(j) Delta2; Delta2' -2*t+mu_F(j)];
        beta_BdG = [-t 0; 0 t];
    
        U=UB;
     
        I_E =  @(E) BTK_IE(E, V(k), t, zplus, ND, mu_F(j), alpha1, alpha2, beta_BdG, kT, U);
        %I_E(i) =  BTK_IE(E(i), V, t, zplus, ND, mu_F, alpha1, alpha2, beta_BdG, kT, U);
        I1_all(u,k) = integral(I_E, -1*abs(V(k)), 1*abs(V(k)),'AbsTol',1e-8,'ArrayValued',true);
    
      toc
    
      end
      k
    end

dV = V(2)-V(1);
G1_all(u,:) = diff(real(I1_all(u,:)))/dV;

end

% % I1_UB(n,:) = real(I1); 

figure(1)
hold on
for u=1:length(UB_array)
   plot(V./Delta2, real(I1_all(u, :))./Delta2, colors(u), 'linewidth', 2);
   % text(V(end)./Delta2, real(I1_all(u, end))./Delta2, labels{u}, 'FontSize', 20, 'HorizontalAlignment', 'right')
end

% title('I-V Characteristics','fontSize',30,'interpreter','latex')
ylabel('Current $(e\Delta_0/\hbar)$','interpreter','latex','fontsize',30)
xlabel('Bias $(eV/\Delta_0)$','fontSize',30,'interpreter','latex')
legend({'U = 0', 'U = 0.5t', 'U = t', 'U = 2t'}, 'Fontsize', 20)
set(gca,'fontSize',30,'linewidth',2,'fontSize',30)

ax = gca;
ax.XAxis.TickDirection = 'in';
ax.YAxis.TickDirection = 'in';
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.XAxis.TickLength = [0.01, 0.005];
ax.YAxis.TickLength = [0.01, 0.005];
ax.FontSize = 20;
ax.Box = 'on';
ax.GridLineStyle = ':';
grid on
grid minor

ylim([0 8])
yticks(0:1:8)
ax.YAxis.MinorTickValues = 0:0.5:8;

xlim([0 2])
xticks(0:0.5:2)
ax.XAxis.MinorTickValues = 0:0.25:2;
  %%

% % figure(2)
% figure(1)
% hold on
% dV = V(2)-V(1);
% for u=1:length(UB_array)
%    plot(V(1:length(V)-1)./Delta2,G1_all(u,:),colors(u),'linewidth', 2);
%     % text(V(end-1)./Delta2, G1_all(u, end), labels{u}, 'FontSize', 20, 'HorizontalAlignment', 'right')
% end
% 
% % title('Conductance','fontSize',30,'interpreter','latex')
% ylabel('G $(e^2/\hbar)$','interpreter','latex','fontsize',30)
% xlabel('Bias $(eV/\Delta_0)$','fontSize',30,'interpreter','latex')
% legend({'U = 0', 'U = 0.5t', 'U = t', 'U = 2t'}, 'Fontsize', 20)
% set(gca,'fontSize',30,'linewidth',2,'fontSize',30)
% 
% ax = gca;
% ax.XAxis.TickDirection = 'in';
% ax.YAxis.TickDirection = 'in';
% ax.XAxis.MinorTick = 'on';
% ax.YAxis.MinorTick = 'on';
% ax.XAxis.TickLength = [0.01, 0.005];
% ax.YAxis.TickLength = [0.01, 0.005];
% ax.FontSize = 20;
% ax.Box = 'on';
% ax.GridLineStyle = ':';
% grid on
% grid minor
% 
% ylim([0 4.5])
% yticks(0:0.5:4.5)
% ax.YAxis.MinorTickValues = 0:0.25:4.5;
% 
% xlim([0 2])
% xticks(0:0.5:2)
% ax.XAxis.MinorTickValues = 0:0.25:2;