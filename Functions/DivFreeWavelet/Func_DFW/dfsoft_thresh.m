% [wcdf1,wcdf2,wcn] = dfsoft_thresh(wcdf1,wcdf2,wcn,minSize,res,FOV,dfthresh,nthresh)
%
% Soft-thresholding operator on divergence-free wavelet coefficients
% 
% Inputs: 
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     FOV         -   image size ( ie size(vx) )
%     dfthresh    -   threshold for divergence-free components
%     nthresh     -   threshold for non-divergence-free components
%
% Outputs:
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%
%