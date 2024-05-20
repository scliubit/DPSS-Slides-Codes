clear all
close all
OMP_iternum = 80;
x = 1:OMP_iternum;
x_sub = 1:2:OMP_iternum;
% load("NMSE_DPSS_dB_acc.mat")
% N = size(NMSE_DPSS_dB,1);
% NMSE_DPSS_dB_mix = zeros(N,1,OMP_iternum);
% NMSE_DPSS_dB_mix(:,:,1:OMP_iternum) = NMSE_DPSS_dB;
load("NMSE_DPSS_dB_acc_long.mat")
NMSE_DPSS_dB_mix = NMSE_DPSS_dB;

figure
hold on
plot(x_sub,squeeze(NMSE_DPSS_dB_mix(1,:,x_sub)),'-*',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB_mix(2,:,x_sub)),'-s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB_mix(3,:,x_sub)),'-v',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB_mix(4,:,x_sub)),'-o',"LineWidth",1.5,'MarkerSize',10)
% plot(x_sub,squeeze(NMSE_DPSS_dB_mix(5,:,x_sub)),'-d',"LineWidth",1.5,'MarkerSize',10)
% plot(x_sub,squeeze(NMSE_DPSS_dB_mix(6,:,x_sub)),'-+',"LineWidth",1.5,'MarkerSize',10)
% plot(x_sub,squeeze(NMSE_DPSS_dB_mix(7,:,x_sub)),'-+',"LineWidth",1.5,'MarkerSize',10)
grid on
legend('Perfect','\epsilon = 0.05','\epsilon = 0.1','\epsilon = 0.2')%'\epsilon = 0.5','\epsilon = 1.0','\epsilon = 2.0')
xlabel('Number of OMP iterations (I)')
ylabel('NMSE [dB]')
axis([1 OMP_iternum -38 0])
if ismac
    set(gca,'fontsize',14);
end
box on

% addpath('breakxaxis/breakxaxis/')
% breakxaxis([75 190])
% rmpath('breakxaxis/breakxaxis/')
