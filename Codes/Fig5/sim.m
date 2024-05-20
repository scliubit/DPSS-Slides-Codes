clear all
close all
addpath('../Fig3/')
%% Para Setup
fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
N = 256;
N_UE = 4;
d = lambda/2;
aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
bandwidth = 9e6;
L=6;
Rician_factor = 1;%1/(1+10);%~10dB
beta = 1;
NT_coord.x = linspace(-aperture/2,aperture/2,N);
NT_coord.y = zeros(1,N);
% H_downlink = zeros(1,N);
% UE location range
UE_range_x = 3;
UE_range_y_min = 3;
UE_range_y_max = 20;

%% DICT Configuration
[Psi_DFT,Psi_POL,Psi_DPSS,index2coor] = dict_design_red(N_UE, N,aperture,k,UE_range_y_min,UE_range_y_max,UE_range_x,beta);%dict_design(N_UE, N,aperture,k,UE_range_y_min,UE_range_y_max,UE_range_x);
%% SIM start
epsilon = 0;
sampleRatio=[0.5];
sampleNumList=round(sampleRatio*N*N_UE);
loc_err = [0,0.05,0.1,0.2,0.5,1.0,2.0];%meter

OMP_iternum = 80;
implementations = 2000;
[aR_DPSS,~] = dpss(N_UE,k*aperture/4/pi/UE_range_y_max/N*N_UE,N_UE,"trace");
NMSE_DFT = zeros(numel(loc_err),implementations,OMP_iternum);
NMSE_POLAR = zeros(numel(loc_err),implementations,OMP_iternum);
NMSE_DPSS = zeros(numel(loc_err),implementations,OMP_iternum);
tic

for err_i = 1:numel(loc_err)
    err = loc_err(err_i)
    sample_i = 1;
    for imp=1:implementations
        UE.x=2*(rand()-0.5)*UE_range_x;
        UE.y=UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
        theta = rand()*2*pi;

%         err = error_HOLO(imp);
        x_err = (rand())*err/2;
        y_err = sqrt(err^2-x_err^2);
        x_hat = UE.x+(-1)^randi(2)*x_err;
        y_hat = UE.y+(-1)^randi(2)*y_err;

%         estUE.x = UE.x+err*rand()*cos(theta);
%         estUE.y = UE.y+err*rand()*sin(theta);
        apperture_UE = (N_UE-1)*d;
        NR_coord.x = UE.x+linspace(-apperture_UE/2,apperture_UE/2,N_UE);
        NR_coord.y = UE.y+zeros(1,N_UE);
        if L>=2
            theta = linspace(0,pi,L-1);
            R = 15;
            scatter_coord.x = R*cos(theta);
            scatter_coord.y = R*sin(theta);
        else
            scatter_coord.x = [];
            scatter_coord.y = [];
        end
        [H_Downlink,H_LoS,H_NLoS,NLoS_steering] = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
%         [H_Downlink,H_LoS,H_NLoS] = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
        H_DL_vec=reshape(H_Downlink,[N*N_UE,1]);

        sampleNum = sampleNumList(sample_i);
        y = zeros(sampleNum,1);
        Phi = zeros(sampleNum,N*N_UE);
        for m = 1:sampleNum
            F_RF = exp(1j*rand(N,4)*2*pi);
            F_BB = (randn(4,1)+randn(4,1)*1j)/sqrt(2);
            f = F_RF*F_BB;
            f_norm = f/norm(f);
            W_RF = exp(1j*rand(N_UE,1)*2*pi);
            W_BB = 1;
            w = W_RF*W_BB;
            w_norm = w'/norm(w);
            y(m) = w_norm*(H_Downlink*f_norm); % ignore noise in the near-field
            Phi(m,:) = kron(f_norm.',w_norm);
        end

%         y_hat = estUE.y;
%         x_hat = index2coor(index,2);
%         x_hat = estUE.x;
        distance_hat = sqrt(y_hat^2+(NT_coord.x-x_hat).^2);
        Psi_DPSS_c = Psi_DPSS.'.*exp(-1j*k.*distance_hat);
        Psi_DPSS_c_multiple = kron(conj(Psi_DPSS_c),aR_DPSS);
        for i = 1:L-1
            NLoS_steering(:,i) = NLoS_steering(:,i)/norm(NLoS_steering(:,i),'fro');
        end
        [ res_sph,~ ] = Successive_OMP( y,Phi,[NLoS_steering,Psi_DPSS_c_multiple'],epsilon,eye(sampleNum,sampleNum),OMP_iternum);
        rec_h_omp_sph = [NLoS_steering,Psi_DPSS_c_multiple']*res_sph;
        for iter_ = 1:OMP_iternum
            res_iter = norm(rec_h_omp_sph(:,iter_)-H_DL_vec)^2/norm(H_DL_vec)^2;
            if isnan(res_iter)
                NMSE_DPSS(err_i,imp,iter_) = 0.001;
            else
                NMSE_DPSS(err_i,imp,iter_) = res_iter;
            end
        end
    end
    toc
end
%% RES Analysis
% index_ = 3; % mu = 0.25
% NMSE_DPSS(:,:,1) = NMSE_POLAR(:,:,1);
% NMSE_DFT_dB = 10*log10(mean(squeeze(NMSE_DFT),2));
% NMSE_POL_dB = 10*log10(mean(squeeze(NMSE_POLAR),2));
NMSE_DPSS_dB = 10*log10(mean(squeeze(NMSE_DPSS),2));
%%
x = 1:OMP_iternum;
x_sub = 1:2:OMP_iternum;
% errorbar(x,squeeze(NMSE_DFT_dB(index,:,:)),err_DFT,'-s',"LineWidth",1.25,'MarkerSize',8)
% hold on
% errorbar(x,squeeze(NMSE_POL_dB(index,:,:)),err_POL,'-s',"LineWidth",1.25,'MarkerSize',8)
% errorbar(x,squeeze(NMSE_SPH_dB(index,:,:)),err_SPH,'-s',"LineWidth",1.25,'MarkerSize',8)
figure
hold on

plot(x_sub,squeeze(NMSE_DPSS_dB(1,:,x_sub)),'-*',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(2,:,x_sub)),'-s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(3,:,x_sub)),'-v',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(4,:,x_sub)),'-o',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(5,:,x_sub)),'-d',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(6,:,x_sub)),'-+',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DPSS_dB(7,:,x_sub)),'-x',"LineWidth",1.5,'MarkerSize',10)
grid on
legend('Perfect','\epsilon = 0.05','\epsilon = 0.1','\epsilon = 0.2','\epsilon = 0.5','\epsilon = 1.0')
xlabel('Number of OMP iterations (I)')
ylabel('NMSE[dB]')
axis([1 OMP_iternum -35 0])
if ismac
    set(gca,'fontsize',14);
end
%% with error bar
% WONTFIX: Difficult to draw
%%
if implementations>=1000
    save('../results/NMSE_DPSS_dB_acc.mat','NMSE_DPSS_dB')
end


rmpath('../Fig3/')