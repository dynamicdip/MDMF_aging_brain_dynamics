function [sim_fc,bold_signal]=simu_fc(T,S_E,ds,nAreas)
% Input:    T     = Total time in seconds
%           S_E   = gating variable. Size: t x N  ; t = time intervals. N = no. of brain areas. 
%           ds    = down-sampling rate in seconds
%           nAreas= no. of brain areas
% Output:  bold_signal of size: t x N.  sim_FC of size N x N.

B = BOLD(T,S_E(:,1)) ; % B=BOLD activity, bf=Fourier transform, f=frequency range)
BOLDAct = zeros(length(B),nAreas);
BOLDAct(:,1) = B ;
for areaNo = 2:nAreas
    BOLDAct(:,areaNo) = BOLD(T,S_E(:,areaNo),0.001);
end
bold_signal = BOLDAct(1:ds:end,:) ;
sim_fc = corrcoef( bold_signal ) ;
sim_fc=sim_fc.*~eye(nAreas);