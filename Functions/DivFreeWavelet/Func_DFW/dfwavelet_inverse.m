% [vx,vy,vz] = dfwavelet_inverse(wcdf1,wcdf2,wcn,minSize,res,FOV)
%
% Inverse Wavelet Transform using divergence-free wavelets
% 
% Inputs: 
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     FOV         -   image size ( ie size(vx) )
%
% Outputs:
%     vx,vy,vz    -   velocity data
%
%
%