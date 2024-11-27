function rE=phie(x)
% excitatory gating variables used in MDMF model
%Input  : excitatory currents, x in nA
%Output : excitatory firing rate, rI in Hz

dE=0.16;
bE=125.;
aE=310.;
y=aE*x-bE;
rE = y./(1-exp(-dE*y));
end