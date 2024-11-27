function [A,BC,CP,E,T,Q,mfc,vfc,kfc,sfc]= return_property_fc_sim(path,thres)

% Inputs: path of the subject wise FCs, e.g.: ".\rok_sc_old.mat"
%        threshold = used to binarize the FCs, e.g., threshold = 0.05;
% Outputs:
% A=assortativity  	[Size: no of subjects x 1 ]
% BC=betweenness   	[Size: no of subjects x 1 ]
% CP= charpath      	[Size: no of subjects x 1 ]
% E=efficiency 	[Size: no of subjects x 1 ]
% T=transitivity 	[Size: no of subjects x 1 ]
% Q=modularity 	[Size: no of subjects x 1 ]


cd(path)
nn  =  dir('*.mat') ;
inx =  1:length(nn) ;
M   =  zeros(length(inx),1) ;
CP  =  zeros(length(inx),1) ;
Q   =  zeros(length(inx),1) ;
GE  =  zeros(length(inx),1) ;
T   =  zeros(length(inx),1) ;
A   =  zeros(length(inx),1) ;
E   =  zeros(length(inx),1) ;
BC  =  zeros(length(inx),1) ;
mfc =  zeros(length(inx),1) ;
vfc =  zeros(length(inx),1) ;
kfc =  zeros(length(inx),1) ;
sfc =  zeros(length(inx),1) ;
i   =  0 ;
for co  = 1:length(inx)
    i   = i+1 ;
    str = nn(inx(co)).name;
    load(str)  ;
    FC=dat(:,:,2);
    mfc(i)=mean(FC(:));
    vfc(i)=var(FC(:));
    kfc(i)=kurtosis(FC(:));
    sfc(i)=skewness(FC(:));
    %%%%----------------prepare FC----------------------------------------
%     FC=abs(f1.FC);
    nAreas = size( FC, 1 ) ; FC( 1 : nAreas+1 : nAreas*nAreas ) = 0 ;
    bU=abs(FC);
    clear FC
    % % % % ---------------binary FC Network-------------------
    bU( bU < thres ) = 0 ; bU( bU >= thres ) = 1 ;
    % %     % % % % ---------------FC analysis-------------------
    M(i)             = mean(modularity_und( bU  )) ;
    [~,Q(i)]         = modularity_und( bU  ) ;
    [CP(i), GE(i) ]  = charpath( bU  )  ;
    E(i)             = efficiency_bin(bU,0);
    T(i)             = transitivity_bu(bU);
    A(i)             = assortativity_bin(bU,0);
    b                = betweenness_bin(bU);
    BC(i)            = mean(b/max(b));
    %
    %     M(i)              = mean(modularity_und( bU  )) ;
    %     [~,Q(i)]          = modularity_und( bU  ) ;
    %     [CP(i), GE(i) ]   = charpath( bU  )  ;
    %     E(i)              = efficiency_wei(bU,0);
    %     T(i)              = transitivity_wu(bU);
    %     A(i)              = assortativity_wei(bU,0);
    %     b                 = betweenness_wei(bU);
    %     BC(i)             = mean(b/max(b));
end

