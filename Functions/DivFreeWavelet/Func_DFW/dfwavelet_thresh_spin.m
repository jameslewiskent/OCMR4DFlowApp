% [vx,vy,vz] = dfwavelet_thresh_spin(vx,vy,vz,minSize,res,dfthresh,nthresh,spins,isRand)
%
% Divergence-free wavelet denoising using specified thresholds with cycle
% spinning
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     dfthresh    -   threshold for divergence-free components
%     nthresh     -   threshold for non-divergence-free components
%     spins       -   log2(number of cycle spinning),(spins=3 equals 8 shifts)
%     isRand      -   specify whether the cycle spinning is random
%
% Outputs:
%     vx,vy,vz    -   denoised velocity data
%
%