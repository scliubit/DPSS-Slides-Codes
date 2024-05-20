clear all
close all
addpath('../no_redundant/')
%% parameters
fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
N = 256;
N_UE = 16;
d = lambda/2;
aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
bandwidth = 9e6;
L=100;
Rician_factor = 0;%1/(1+10);%~10dB

NT_coord.x = linspace(-aperture/2,aperture/2,N);
NT_coord.y = zeros(1,N);
UE_range_x = 2;
UE_range_y_min = 2;
UE_range_y_max = 20;
x_lim = (N-1)*d;
x_range = linspace(-x_lim/2,x_lim/2,N);


UE.x=0;%2*(rand()-0.5)*UE_range_x;
UE.y=2;%UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
apperture_UE = (N_UE-1)*d;
NR_coord.x = UE.x+linspace(-apperture_UE/2,apperture_UE/2,N_UE);
NR_coord.y = UE.y+zeros(1,N_UE);
if L>=2
    theta = linspace(-pi,pi,L-1);
    R = 25;
    scatter_coord.x = R*cos(theta);
    scatter_coord.y = R*sin(theta);
else
    scatter_coord.x = [];
    scatter_coord.y = [];
end
[H_Downlink,H_LoS,H_NLoS] = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
H_downlink=H_Downlink(1,:);
figure

%% DFT
Dict01 = dftmtx(N)/sqrt(N);
%% POL
UE_range_z_max = 5;
UE_range_z_min = 1;
UE_range_x = 2;
r_inv = linspace(1/UE_range_z_max,1/UE_range_z_min,round(sqrt(N)));%0:0.01:1;
r = 1./r_inv;
x = linspace(-UE_range_x,UE_range_x,round(sqrt(N)));
DICT_Layer1 = zeros(N,numel(r_inv));
cnt=0;
index2coor=zeros(numel(r_inv),2);
for xx =1:numel(x)
    for rr = 1:numel(r)
        cnt=cnt+1;
        %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
        DICT_Layer1(:,cnt)=exp(1j*k.*sqrt(r(rr)^2+(x(xx)-x_range.').^2))/sqrt(N);%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
        index2coor(cnt,1) = r(rr);
        index2coor(cnt,2) = x(xx);
    end
end
DICT_Layer1=DICT_Layer1(:,1:N);
%% prolate spheroidal (Slepian) functions
seq_length = N; 
time_halfbandwidth = 4;%k*aperture/4/pi/UE_range_z_max;
% num_seq = 4;
[dps_seq,lambda] = dpss(N,time_halfbandwidth,N,"trace");

H_downlink = zeros(1,numel(x_range));
recv_x = 1;
z = 2;

distance = sqrt(z^2+(x_range-recv_x).^2);
recv_em = zeros(1,1);
dps_seq_c = dps_seq'.*exp(-1j*k.*distance);%exp(-1j*k.*(x_range).^2/z);
% dps_seq_c = dps_seq'.*exp(-1j*k.*(x_range).^2/z/2);

for i=1:numel(x_range)
    recv_em = recv_em+1/distance(i)*exp(-1j*k*distance(i));
    H_downlink(i) = 1/distance(i)*exp(-1j*k*distance(i));
end

subplot(311)
hold on
stem(abs(Dict01*H_downlink'),'LineWidth',1.25,'MarkerSize',10,LineStyle='-')
stem(abs(DICT_Layer1'*H_downlink'),'LineWidth',1.25,'MarkerSize',10,LineStyle='--',Marker='v')
stem(abs(dps_seq_c*H_downlink'),'LineWidth',1.25,'MarkerSize',10,LineStyle='--',Marker='square')
grid on
axis([0 N 0 5.2])
xlabel({'Dictionary Index';'(b) Mutual Correlation'},'Interpreter','latex')
ylabel('Coefficient','Interpreter','latex')
legend({'DFT Dictionary','Spherical Dictionary','Proposed'},'Interpreter','latex')
title('(a) Sparse Representations','Interpreter','latex')
if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','times new roman');
end
box on

%% ortho

subplot(323)

surf(10*log10(1+fliplr(abs(DICT_Layer1'*DICT_Layer1)))/3-0.001)
cmp = colormap('gray');
cmp = flipud(cmp);
colormap(cmp);
colorbar('XTickLabel',{'0','1'},'XTick', 0:1)
shading interp
% ax.YAxis.Visible = 'off';
% ax.XAxis.Visible = 'off';
view([0 0 1])
axis('square')
axis('tight')
xticks([128,256])
yticks([128,256])
% YTick(['50','100','150'])
xlabel('Index','Interpreter','latex')
ylabel('Index','Interpreter','latex')
title('Spherical','Interpreter','latex')
if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','times new roman');
end
box on


subplot(324)

surf(10*log10(1+fliplr(abs(dps_seq_c'*dps_seq_c)))/3)
cmp = colormap('gray');
cmp = flipud(cmp);
colormap(cmp);
colorbar('XTickLabel',{'0','1'},'XTick', 0:1)
shading interp
view([0 0 1])
axis('square')
axis('tight')
xticks([128,256])
yticks([128,256])
% YTick(['50','100','150'])
xlabel('Index','Interpreter','latex')
ylabel('Index','Interpreter','latex')
title('Proposed','Interpreter','latex')

if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','times new roman');
end
box on


%% sinc
% addpath('../no_redundant/')
RT = H_Downlink'*H_Downlink;
RR = H_Downlink*H_Downlink';
% figure;surf(real(RT));shading interp
x=0:d:aperture;
y = 2*UE.y*sin(k*apperture_UE/2/UE.y.*x)./x/k/apperture_UE;
subplot(313)
hold on
plot(1:N,abs(y),LineWidth=1.5)
sep = 4;
x_down = 1:sep:N;
plot(x_down,abs(RT(1,1:sep:N))./abs(RT(1,1)),'o',LineWidth=1.5,MarkerSize=8)
axis([0 N 0 1])
legend({'Sinc Function','Auto-Correlation'},'Interpreter','latex')
xlabel('Antenna Index','Interpreter','latex')
ylabel({'Normalized';'Auto-correlation'},'Interpreter','latex')
title('(c) Approximation Error','Interpreter','latex')
grid on
% figure;surf(real(RR));shading interp
% rmpath('../no_redundant/')

if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','times new roman');
end
box on

set(gcf,'renderer','Painters')