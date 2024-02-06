% [vx,vy,vz] = dfwavelet_thresh_SURE_MAD(vx,vy,vz,minSize,res)
%
% Divergence-free wavelet denoising using Sure-based threhsold and MAD
% sigma estimation
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%
% Outputs:
%     vx,vy,vz    -   denoised velocity data
%
%