% [wcdf1,wcdf2,wcn,numLevels,wc_sizes] = dfwavelet_forward(vx,vy,vz,minSize,res)
%
% Forward Wavelet Transform using divergence-free wavelets
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%
% Outputs:
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%     numLevels   -   Number of wavelet decomposition levels
%     wc_sizes    -   wavelet subband sizes 
%                     (The first three numbers correspond to the scaling
%                     subband, the second group corresponds to the smallest
%                     detail subband and it gets bigger and bigger
%
%


