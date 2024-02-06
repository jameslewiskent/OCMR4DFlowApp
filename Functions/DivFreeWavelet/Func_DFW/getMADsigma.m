% sigma = getMADsigma(wcn,minSize,res,FOV,percentZero)
%
% Estimate noise standard deviation using the Median Absolute Deviation on
% the non-divergence-free wavelet subband
% 
% Inputs: 
%     wcn         -   Non-divergence-free wavelet coefficients
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     FOV         -   image size ( ie size(vx) )
%     percentZero -   percentage of masked out elements in velocity data
%
% Outputs:
%     sigma       -   noise standard deviation
%
%