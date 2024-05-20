function [H_Downlink,H_LoS,H_NLoS,NLoS_steering] = GenChannelDL(NR,NT,L,NR_coord,NT_coord,scatter_coord,kappa,Rician_factor)
H_Downlink = zeros(NR,NT);
H_LoS = zeros(NR,NT);
H_NLoS = zeros(NR,NT);
assert(L>=1); % number of multipaths
assert(numel(NR_coord.x)==NR)
assert(numel(NT_coord.x)==NT)
assert(numel(scatter_coord.x)==L-1)
NLoS_steering = zeros(NR*NT,L-1);
for ell = 1:L
    if ell == 1
        % generate LoS path
        for n = 1:NR
            distance = sqrt((NT_coord.x-NR_coord.x(n)).^2+(NT_coord.y-NR_coord.y(n)).^2);
            H_LoS(n,:) = 1./distance.*exp(-1j*kappa.*distance);
        end
%         pwr = norm(H_LoS,'fro');
    else
        % break
        distance1 = sqrt((NT_coord.x-scatter_coord.x(ell-1)).^2+(NT_coord.y-scatter_coord.y(ell-1)).^2);
        distance2 = sqrt((scatter_coord.x(ell-1)-NR_coord.x).^2+(scatter_coord.y(ell-1)-NR_coord.y).^2)';
        NLoS_steering(:,ell-1) = kron((1./distance1.*exp(-1j*kappa.*distance1)).',(1./distance2.*exp(-1j*kappa.*distance2)));
%         H_NLoS = H_NLoS + (randn()+1j*randn())/sqrt(2)*pwr*sqrt(Rician_factor/(L-1))*kron(1./distance2.*exp(-1j*kappa.*distance2),1./distance1.*exp(1j*kappa.*distance1));
%         H_NLoS = H_NLoS + (randn()+1j*randn())/sqrt(2)*sqrt(Rician_factor)*((1./distance2.*exp(-1j*kappa.*distance2))*(1./distance1.*exp(1j*kappa.*distance1)));
%         H_NLoS = H_NLoS + (randn()+1j*randn())/sqrt(2)*(1/L-1)*(1./distance2.*exp(-1j*kappa.*distance2))*(1./distance1.*exp(-1j*kappa.*distance1));
        H_NLoS = H_NLoS + sqrt(1/L-1)*(1./distance2.*exp(-1j*kappa.*distance2))*(1./distance1.*exp(-1j*kappa.*distance1));
    end
end
% pwr_n = norm(H_NLoS,'fro');
% pwr/pwr_n
% H_NLoS = sqrt(Rician_factor)*(randn(NR,NT)+randn(NR,NT)*1j)/sqrt(2)/sqrt(NR*NT);
H_Downlink = H_LoS+sqrt(Rician_factor)*H_NLoS;
end

