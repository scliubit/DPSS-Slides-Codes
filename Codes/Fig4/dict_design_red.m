function [Psi_DFT,Psi_POL,Psi_DPSS,index2coor] = dict_design_red(N_UE,N,aperture,k,z_min,z_max,x_max,beta)

x_range = linspace(-aperture/2,aperture/2,N);
aperture_UE = aperture/N*N_UE;
% DFT DICT
if N_UE>1
    aR_DFT = dftmtx(N_UE*beta)/sqrt(N_UE*beta);
    aR_DFT = aR_DFT(1:N_UE,:);
    aT_DFT = dftmtx(N*beta)/sqrt(N*beta);
    aT_DFT = aT_DFT(1:N,:);
    Psi_DFT = kron(conj(aT_DFT),aR_DFT);
else
    Psi_DFT = dftmtx(N*beta)/sqrt(N*beta);
    Psi_DFT = Psi_DFT(1:N,:);
end
% POL DICT
if N_UE>1
    r_inv = linspace(1/z_max,1/z_min,floor(sqrt(N*beta)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,ceil(sqrt(N*beta)));
    aT_POL = zeros(numel(x_range),numel(r_inv)*numel(x));
    cnt=0;
    index2coor=zeros(numel(r_inv)*numel(x),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aT_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = r(rr);
            index2coor(cnt,2) = x(xx);
        end
    end
    index2coor = kron(index2coor,[1;1;1;1]);
    % ==============================
    x_range_R = linspace(-aperture_UE/2,aperture_UE/2,N_UE);
    r_inv = linspace(1/z_max,1/z_min,floor(sqrt(N_UE)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,ceil(sqrt(N_UE)));
    aR_POL = zeros(numel(x_range_R),numel(r_inv)*numel(x));
    cnt=0;
    index2coor2=zeros(numel(r_inv),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aR_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range_R.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor2(cnt,1) = r(rr);
            index2coor2(cnt,2) = x(xx);
        end
    end
    Psi_POL = kron(aT_POL,aR_POL);
else
    r_inv = linspace(1/z_max,1/z_min,floor(sqrt(N*beta)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,ceil(sqrt(N*beta)));
    Psi_POL = zeros(numel(x_range),numel(r_inv)*numel(x));
    cnt=0;
    index2coor=zeros(numel(r_inv),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            Psi_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = r(rr);
            index2coor(cnt,2) = x(xx);
        end
    end
end
% Prolate DICT


if N_UE>1
    time_halfbandwidth = 0;%k*aperture/4/pi/z_max;
    time_halfbandwidth_R = 1;%k*aperture_UE/4/pi/z_max;
    % num_seq = 4;
    [aT_DPSS,~] = dpss(N,time_halfbandwidth,N,"trace");
    [aR_DPSS,~] = dpss(N_UE,time_halfbandwidth_R,N_UE,"trace");
    Psi_DPSS = aT_DPSS;%kron(conj(aT_DPSS),aR_DPSS);
else
    time_halfbandwidth = 4;%k*aperture/4/pi/z_max;
    % num_seq = 4;
    [Psi_DPSS,~] = dpss(N,time_halfbandwidth,N,"trace");
end

end

