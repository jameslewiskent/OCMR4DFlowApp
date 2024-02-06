% [vx,vy,vz] = dfwavelet_thresh_SURE(vx,vy,vz,minSize,res,sigma)
%
% Divergence-free wavelet denoising using Sure-based threhsold
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     sigma       -   noise standard deviation
%
% Outputs:
%     vx,vy,vz    -   denoised velocity data
%
%
%