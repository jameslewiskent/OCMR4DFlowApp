%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [wcdf1,wcdf2,wcn] = dfSUREshrink(wcdf1,wcdf2,wcn,minSize,res,FOV,percentZero,sigma)
%
% SUREshrink operator on divergence-free wavelet coefficients
% 
% Inputs: 
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     FOV         -   image size ( ie size(vx) )
%     percentZero -   percentage of masked out elements in velocity data
%     sigma       -   noise standard deviation
%
% Outputs:
%     wcdf1,wcdf2 -   Two divergence-free wavelet coefficients
%     wcn         -   Non-divergence-free wavelet coefficients
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vx,vy,vz] = dfwavelet_thresh(vx,vy,vz,minSize,res,dfthresh,nthresh)
%
% Divergence-free wavelet denoising using specified thresholds
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     dfthresh    -   threshold for divergence-free components
%     nthresh     -   threshold for non-divergence-free components
%
% Outputs:
%     vx,vy,vz    -   denoised velocity data
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vx,vy,vz] = dfwavelet_thresh_SURE_spin(vx,vy,vz,minSize,res,sigma,spins,isRand)
%
% Divergence-free wavelet denoising using Sure-based threhsold and cycle
% spinning
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vx,vy,vz] = dfwavelet_thresh_SURE_MAD_spin(vx,vy,vz,minSize,res,spins,isRand)
%
% Divergence-free wavelet denoising using Sure-based threhsold, MAD and
% cycle spinning
% 
% Inputs: 
%     vx,vy,vz    -   velocity data
%     minSize     -   minimum size for the wavelet scaling subband
%     res         -   resolution (a 1x3 array)
%     spins       -   log2(number of cycle spinning),(spins=3 equals 8 shifts)
%     isRand      -   specify whether the cycle spinning is random
%
% Outputs:
%     vx,vy,vz    -   denoised velocity data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%
% (c) Frank Ong 2013




