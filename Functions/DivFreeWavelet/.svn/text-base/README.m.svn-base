% This is a collection of Matlab functions and demos to reproduce some of
% the results that are described in the paper:
% F Ong, M Uecker, U Tariq, A Hsiao, MT Alley, SS Vasanawala and M Lustig
%"Robust 4D Flow Denoising using Divergence-free Wavelet Transform"
%
%
% You are encouraged to modify/distribute this code. However, please 
% acknowledge this code and cite the papers appropriately. For any 
% questions about the code, please contact me at frankong@berkeley.edu. 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Setup:
%
% To compile the mex files for your architecture, execute the make.m script
%
% Run setpath.m to add all the paths before running any demos
%
% After running make.m and setpath.m, try demo_couette.m for a quick and
% simple demonstration. For other demonstrations, please download
% the DivFreeWavelet_data zip file.
%
% The demos are recommended to run using the Matlab publish tool, which
% produces beautiful documents
%
% Type "help function" to see documentation for the particular function
% See DFW_help.m for a collection of divergence-free wavelet function
% documentation
%
%
% Demos:
%
%       demo_couette: 
%               Compares divergence-free wavelet denoising with existing 
%               methods on a simulated Couette flow. This demo is fast and
%               is recommended to try first before other experiments. This
%               test is not in the MRM paper.
%
%       demo_cfd: 
%               Compares divergence-free wavelet denoising with existing 
%               methods on a flow data generated in OpenFOAM. 
%
%       demo_phantom: 
%               Compares divergence-free wavelet denoising with existing 
%               methods on a flow phantom data with a correct segmentation
%
%       demo_phantom_seg: 
%               Compares divergence-free wavelet denoising with existing 
%               methods on a flow phantom data with a incorrect segmentation
%
%       demo_phantom_subsample: 
%               Compares divergence-free wavelet denoising with existing 
%               methods on an undersampled flow phantom data reconstructed
%               with ESPIRiT/Iterative SENSE with no regularization
%
%
%
% General Comments about the code:
%
%       All velocity data are oriented like a matrix. That is vx(x,y,z)
%       corresponds to vx(i,j,k) directly
%       Positive x direction corresponds to "down", positive y direction 
%       corresponds to "right" and positive z direction corresponds to "out".
%
%       All plotting functions have a slider for slice selection. But if
%       the slice is specified as an input, all controls will not appear,
%       which makes it better to export the figures.
%      
%       Mex functions for DFW doesn't work well with OpenMP and CUDA (in Mac)
%       so they are compiled without them by default.
%
%       Within ./Func_DFW/src, the C and CUDA files can be used
%       separately from the mex interfaces. The default MakeFile enables 
%       OpenMP but not CUDA. 
%       To compile CUDA, set CUDA=1 and change cuda.top to the correct CUDA
%       library path. Run ./dfwavelet to see if the files are
%       correctly compiled.
%
%
% 
% For any questions, comments and contributions, please contact
% Frank Ong (frankong@berkeley.edu)
% 
% 

