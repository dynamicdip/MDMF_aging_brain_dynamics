function [ FCD_se, R_se, std_R_se, FCD_re, R_re, std_R_re ] = MDMF_fr( G, T_Glu, T_Gaba, SC, fcEmp )

% Inputs:G     = global coupling strength % e.g., G=0.75
%        T_Glu  = the glutamate level    % e.g., T_Glu = 11.3
%        T_Gaba = the GABA level    % e.g., T_Gaba = 3.35
%        SC     = structural connectivity matrix drawn from DTI of size: (area  x area)
%        FC     = empirical functional connectivity matrix derived from empirical BOLD signals. Size: (area  x area)
%  Outputs: FCD = simulated functional connectivity matrix from simulated BOLD signals. Size: (area  x area)
%           R_se = synchrony index , std_R_se = metastability index.

nAreas  = size(SC,1);                 % number of brain areas determined based on structural connectivity matrix, Craddock parcellation is used
sigma   = 0.001;               % noise intensity
dt      = 0.01;		       % time interval for simulation
simTime = 5*60*1000;           % total simulation time
tSpan   = 0:dt:simTime;
tr      = 2*60*1000;           % transient time. Discard this amount of simulation.
SE = 0.01*rand(nAreas,1) ;     % initial values
SI = 0.01*rand(nAreas,1) ;     % initial values
J  = 0.1*rand(nAreas,1) ;      % initial values
% %% Model parameters
JNMDA = 0.15;       % scaling factor to get low spontaneous activity for isolated node
I0 = 0.382;         % input current to each node
WE = 1.0;           % scaling of external input current to excitatory population
WI = 0.7;           % scaling of external input current to inhibitory population
wplus = 1.4;        % weight for recurrent self-excitation in each excitatory population
gamma = 1.0;        % learning rate
alphaE = 0.072;     % forward rate constant for NMDA gating
betaE = 0.0066;     % backward rate constant for NMDA gating
alphaI = 0.53;      % forward rate constant for GABA gating
betaI = 0.18;       % backward rate constant for GABA gating
rho = 3;            % target-firing rate of the excitatory population is maintained at the 3 Hz
% %% Initializing Variables
S_E = zeros(simTime,nAreas);        % initialization of excitatory synaptic gating variable
S_I = zeros(simTime,nAreas);        % initialization of inhibitory synaptic gating variable
J_IE = zeros(simTime,nAreas);
r_E = zeros(simTime,nAreas);
simOut.maxfr_E = zeros(nAreas,1);   % initialization of mean firing rate for excitatory population of each brain area
simOut.maxfr_I = zeros(nAreas,1);   % initialization of mean firing rate for inhibitory population of each brain area
% %% --------------------------------
nn = 1 ; j = 0 ;
for i  = 2:1:length(tSpan)
    IE = WE*I0 + wplus*JNMDA*SE + G*JNMDA*SC*SE-J.*SI;          % input current to the excitatory population
    II = WI*I0 + JNMDA*SE - SI;                                 % input current to the inhibitory population
    rE = phie(IE);                                              % firing rate of excitatory population in brain area
    rI = phii(II);                                              % firing rate of inhibitory population in brain area
    SEdot = - betaE*SE+(alphaE./1000.0)*T_Glu*(1-SE).*rE ;      % kinetic equation represents glutamate binding with recepor (NMDA)
    nun   = randn(nAreas,1);
    SE    = SE + dt*SEdot + sqrt(dt)*sigma*nun ;                % noise term added (uncorrected standard Gaussian noise with the noise amplitude (sigma) in each brain area)
    SE(SE>1) = 1;    SE(SE<0) = 0;
    nug   = randn(nAreas,1);
    SIdot = -betaI*SI+(alphaI./1000.0)*T_Gaba*(1-SI).*rI;       % kinetic equation represents GABA binding with inhibitory receptor
    SI    = SI + dt*SIdot + sqrt(dt)*sigma*nug;                 % noise term, uncorrected standard Gaussian noise with the noise amplitude (sigma) in each brain area)
    SI(SI>1) = 1;    SI(SI<0) = 0;
    Jdot  = gamma*(rI./1000.0).*((rE-rho)./1000.0);              % J, synaptic weight and its dynamics clamps firing rate of the excitatory population
    J     = J + dt*Jdot ;
    j=j+1;
    if j==1/dt        
        S_E(nn,:)  = SE' ;
        S_I(nn,:)  = SI' ;        
        J_IE(nn,:) = J'  ;
        r_E(nn,:) =  rE' ;
        if(nn>tr)
            simOut.maxfr_E = max(simOut.maxfr_E,rE);
            simOut.maxfr_I = max(simOut.maxfr_I,rI);
        end
        nn=nn+1;        j=0;
    end    
end
% disp('simulation completed and synaptic activity computed')
nn = nn-1;
% %% Parameters controlling numerics
dtt   = 1e-3;         % sampling rate of simulated neuronal activity (seconds)
ds    = 100;          % BOLD downsampling rate
T = nn*dtt - 2*60 ; % Total time in seconds ignoring the first 2*60 seconds to allow for initial transients
se=S_E(tr+1:end,:);
re=r_E(tr+1:end,:);
sim_fc_se=simu_fc(T,se,ds,nAreas);
sim_fc_re=simu_fc(T,re,ds,nAreas);
FCD_se = fc_distance( sim_fc_se, fcEmp ) ;
FCD_re = fc_distance( sim_fc_re, fcEmp ) ;
% %% ----------------------------------------------
% HighPass = 0.06 ; LowPass = 0.03 ; fs = 1/2 ;
% BP_simulatedMtrx = BPmeta( BOLDActDS, HighPass, LowPass, fs ) ;
% [ std_R, R ] = MetaStab( BP_simulatedMtrx ) ;
% %% ----------------------------------------------
[R_se,std_R_se]=meta_stability(se);
[R_re,std_R_re]=meta_stability(re);