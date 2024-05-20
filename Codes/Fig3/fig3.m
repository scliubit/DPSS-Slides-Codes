clear all
close all
load('../results/NMSE_DFT_dB.mat')
load('../results/NMSE_POL_dB.mat')
load('../results/NMSE_DPSS_dB.mat')
OMP_iternum = 45;
x = 1:OMP_iternum;
x_sub = 1:2:OMP_iternum;
% errorbar(x,squeeze(NMSE_DFT_dB(index,:,:)),err_DFT,'-s',"LineWidth",1.25,'MarkerSize',8)
% hold on
% errorbar(x,squeeze(NMSE_POL_dB(index,:,:)),err_POL,'-s',"LineWidth",1.25,'MarkerSize',8)
% errorbar(x,squeeze(NMSE_SPH_dB(index,:,:)),err_SPH,'-s',"LineWidth",1.25,'MarkerSize',8)
figure
plot(x_sub,squeeze(NMSE_DFT_dB(1,:,x_sub)),'-o',"LineWidth",1.5,'MarkerSize',10)
hold on
plot(x_sub,squeeze(NMSE_DFT_dB(2,:,x_sub)),'--o',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DFT_dB(3,:,x_sub)),':o',"LineWidth",1.5,'MarkerSize',10)


plot(x_sub,squeeze(NMSE_POL_dB(1,:,x_sub)),'-v',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_POL_dB(2,:,x_sub)),'--v',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_POL_dB(3,:,x_sub)),':v',"LineWidth",1.5,'MarkerSize',10)

plot(x_sub,smooth(squeeze(NMSE_DPSS_dB(1,:,x_sub)),7),'-s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,smooth(squeeze(NMSE_DPSS_dB(2,:,x_sub)),7),'--s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,smooth(squeeze(NMSE_DPSS_dB(3,:,x_sub)),7),':s',"LineWidth",1.5,'MarkerSize',10)

grid on
legend({'DFT Dictionary, $\mu=0.25$','DFT Dictionary, $\mu=0.4$','DFT Dictionary, $\mu=0.6$','Spherical Dictionary, $\mu=0.25$','Spherical Dictionary, $\mu=0.4$','Spherical Dictionary, $\mu=0.6$','Proposed, $\mu=0.25$','Proposed, $\mu=0.4$','Proposed, $\mu=0.6$'},'Interpreter','Latex','Location','best')
xlabel('Number of OMP Iterations ($\it I$)','Interpreter','Latex')
ylabel('NMSE [dB]','Interpreter','Latex')
axis([1 OMP_iternum-1 -25 0])
% if ismac
%     set(gca,'fontsize',14);
% end
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on