function [zR,zstd_R]=meta_stability(se)
%[R,std_R]=meta_stability(bold_signal)
%Input: BOLD time series,  t x N.
%Output: R=synchronization index and std_R= metastability index. Outputs are two scaler values.


zsse = (se-mean(se))./std(se) ;
phs  = angle(hilbert(zsse)) ;
phs  = exp(1i*phs);          % complex coordinate
mphs_Rt = abs(mean(phs,2)) ;% average over nodes
zstd_R  = std(mphs_Rt) ;      % standard deviation, metastability
zR      = mean(mphs_Rt) ;