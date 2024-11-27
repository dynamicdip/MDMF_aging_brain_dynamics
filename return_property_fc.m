function [A,BC,CP,E,T,Q]= return_property_fc(path,thres)
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
i   =  0 ;
for co  = 1:length(inx)
    i   = i+1 ;
    str = nn(inx(co)).name;
    f1=load(str)  ;

    %%%%----------------prepare FC----------------------------------------
    FC=abs(f1.FC);  nAreas = size( FC, 1 ) ; FC( 1 : nAreas+1 : nAreas*nAreas ) = 0 ;
    bU=FC;
    clear FC
    % % % % ---------------binary FC Network-------------------
    th = thres;    % th = 0.01 ;
    bU( bU < th ) = 0 ; bU( bU >= th ) = 1 ;
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

