function rI=phii(x)
% inhibitory gating variables used in MDMF model
% Input  : inhibitory currents, x in nA
% Output : inhibitory firing rate, rI in Hz
dI=0.087;
bI=177.;
aI=615.;
y=aI*x-bI;
rI = y./(1-exp(-dI*y));
end