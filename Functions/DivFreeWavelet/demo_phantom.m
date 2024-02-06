%% Denoise Flow Phantom Data
%
% The following compares divergence-free wavelet denoising with existing
% methods on a flow phantom with correct segmentation.
%
% see demo_phantom_seg.m for performances with an incorrect segmentation

%% Clear all and check path
close all
clc
clear

if ~exist('dfwavelet_thresh','file')
    error('Cannot find dfwavelet functions. run setPath!');
end

%% Load phantom data and set parameters
load 4dFlow_007


% Eddy Current Correction 
calib = [127-100,128+100;127-40,128+40;32-20,32+20];
[vxEddy,vyEddy,vzEddy] = eddyCorr(vx,vy,vz,calib);
vx = vx-vxEddy;
vy = vy-vyEddy;
vz = vz-vzEddy;


% Crop to pipe for faster performance
% Comment the following to denoise the entire flow data
imMag = imMag(20:250,1:64,17:49);
vx = vx(20:250,1:64,17:49);
vy = vy(20:250,1:64,17:49);
vz = vz(20:250,1:64,17:49);
vxEddy = vxEddy(20:250,1:64,17:49);
vyEddy = vyEddy(20:250,1:64,17:49);
vzEddy = vzEddy(20:250,1:64,17:49);


% Set parameters
imThresh = 7000;                % Image threshold for segmentation
imMask = (imMag>imThresh)*1;    % Segmentation mask
vMag = getVelMag(vx,vy,vz);     % Velocity magnitude
vMax = 1.0127e+03;              % Maximum speed
FOV = size(vMag);               % Field of view
N = FOV(3);                     % Number of slices (for plotting)
ph0 = zeros(FOV);               % Reference phase in phase contrast


% Plot
figure,imshow_flow(imMask,vx,vy,vz,vMax,[1,2,3]),title('Original Flow Field','FontSize',14);
figure,imshow_flowmag(imMask,vx,vy,vz,vMax,[1,2,3]),title('Original Velocity Magnitude','FontSize',14);
pause(1);

%% Get Noisy Data
load 4dFlow_007_noisy


% Crop to pipe for faster performance
% Comment the following to denoise the entire flow data
vxN = vxN(20:250,1:64,17:49);
vyN = vyN(20:250,1:64,17:49);
vzN = vzN(20:250,1:64,17:49);


% Eddy Current Correction 
vxN = vxN-vxEddy;
vyN = vyN-vyEddy;
vzN = vzN-vzEddy;


% Mask flow data with good segmentation
[vxN,vyN,vzN] = maskIM(imMask,vxN,vyN,vzN);


% Plot
figure,imshow_flow(imMask,vxN,vyN,vzN,vMax,[1,2,3]),title('Noisy Flow Field','FontSize',14);
figure,imshow_flowmag(imMask,vxN,vyN,vzN,vMax,[1,2,3]),title('Noisy Velocity Magnitude','FontSize',14);


% Calculate errors
disp('Noisy Flow Error')
[vNRMSE_Noise,vMagErr_Noise,angErr_Noise] = calcVelError(imMask,vx,vy,vz,vxN,vyN,vzN);
PVNR = 20*log10(1/vNRMSE_Noise);
fprintf('PVNR: \t\t\t%.2fdB\nNRMSE: \t\t\t%f\nvMag Error: \t\t%f\nAbsolute Angle Error: \t%f\n\n',PVNR,vNRMSE_Noise,vMagErr_Noise,angErr_Noise);

pause(1)

%% DivFree Wavelet with SureShrink and MAD sigma estimation
% Here, we use Median Absolute Deviation to estimate noise std
% and then use SureShrink to find the optimal threshold that minimizes MSE


minSize = 8*ones(1,3); % Smallest wavelet level size

% Denoise
[vxDFWsm,vyDFWsm,vzDFWsm] = dfwavelet_thresh_SURE_MAD(vxN,vyN,vzN,minSize,res);


% Plot
figure,imshow_flow(imMask,vxDFWsm,vyDFWsm,vzDFWsm,vMax,[1,2,3])
title('Div Free Wavelet w/ SureShrink (Flow Field)','FontSize',14)
figure,imshow_flowmag(imMask,vxDFWsm,vyDFWsm,vzDFWsm,vMax,[1,2,3])
title('Div Free Wavelet w/ SureShrink (Vel Mag)','FontSize',14)


% Calculate errors
disp('DivFree Wavelet w/ SureShrink and MAD')
[vNRMSE_DFWsm,vMagErr_DFWsm,angErr_DFWsm] = calcVelError(imMask,vx,vy,vz,vxDFWsm,vyDFWsm,vzDFWsm);
fprintf('NRMSE: \t\t\t%f\nvMag Error: \t\t%f\nAbsolute Angle Error: \t%f\n\n',vNRMSE_DFWsm,vMagErr_DFWsm,angErr_DFWsm);
pause(1);

%% DivFree Wavelet with SureShrink, MAD and random cycle spinning
% To remove the blocking artifacts, we do partial cycle spinning
% Here we do 2^3=8 random shifts


spins = 2;              % Number of cycle spinning per dimension
isRandShift = 1;        % Use random shift
minSize = 8*ones(1,3);  % Smallest wavelet level size


% Denoise
[vxDFWsms,vyDFWsms,vzDFWsms] = dfwavelet_thresh_SURE_MAD_spin(vxN,vyN,vzN,minSize,res,spins,isRandShift);


% Plot
figure,imshow_flow(imMask,vxDFWsms,vyDFWsms,vzDFWsms,vMax,[1,2,3])
title('Div Free Wavelet w/ SureShrink and Partial Cycle Spinning (Flow Field)','FontSize',14)
figure,imshow_flowmag(imMask,vxDFWsms,vyDFWsms,vzDFWsms,vMax,[1,2,3])
title('Div Free Wavelet w/ SureShrink and Partial Cycle Spinning (Vel Mag)','FontSize',14)


% Calculate errors
disp('DivFree Wavelet w/ SureShrink, MAD and Partial Cycle Spinning')
[vNRMSE_DFWsms,vMagErr_DFWsms,angErr_DFWsms] = calcVelError(imMask,vx,vy,vz,vxDFWsms,vyDFWsms,vzDFWsms);
fprintf('NRMSE: \t\t\t%f\nvMag Error: \t\t%f\nAbsolute Angle Error: \t%f\n\n',vNRMSE_DFWsms,vMagErr_DFWsms,angErr_DFWsms);
pause(1);

%% Finite Difference Method
% The following implements finite difference method denoising 
% as described in:
% Song SM, Pelc NJ., et al. JMRI 1993
% Noise reduction in three-dimensional phase-contrast MR velocity measurements.


% Denoise
[vxFDM,vyFDM,vzFDM] = fdmDenoise(vxN,vyN,vzN,res);


% Plot
figure,imshow_flow(imMask,vxFDM,vyFDM,vzFDM,vMax,[1,2,3])
title('Div Free FDM','FontSize',14)
figure,imshow_flowmag(imMask,vxFDM,vyFDM,vzFDM,vMax,[1,2,3])
title('Div Free FDM','FontSize',14)


% Calculate errors
disp('Finite Difference Method Error')
[vNRMSE_FDM,vMagErr_FDM,angErr_FDM] = calcVelError(imMask,vx,vy,vz,vxFDM,vyFDM,vzFDM);
fprintf('NRMSE: \t\t\t%f\nvMag Error: \t\t%f\nAbsolute Angle Error: \t%f\n\n',vNRMSE_FDM,vMagErr_FDM,angErr_FDM);
pause(1);

%% DivFree Radial Basis Function
% The following implements divergence-free radial basis function denoising 
% as described in:
% Busch J, Kozerke S., et al. MRM 2012
% Construction of divergence-free velocity fields from cine 3D phase-contrast flow measurements. 


radius = 4;     % Radius of kernel
nIter = 20;     % Number of iterations for lsqr

% Plot during iterations, if on, does gradient descent instead of lsqr
doplot = 0; 


% Denoise
% Takes a while to run
[vxRBF,vyRBF,vzRBF] = rbfDenoise(vxN,vyN,vzN,imMask,radius,res,nIter,doplot);


% Plot
figure,imshow_flow(imMask,vxRBF,vyRBF,vzRBF,vMax,[1,2,3]),
title('Div Free RBF','FontSize',14)
figure,imshow_flowmag(imMask,vxRBF,vyRBF,vzRBF,vMax,[1,2,3]),
title('Div Free RBF','FontSize',14)


% Calculate errors
disp('Radial Basis Function Error')
[vNRMSE_RBF,vMagErr_RBF,angErr_RBF] = calcVelError(imMask,vx,vy,vz,vxRBF,vyRBF,vzRBF);
fprintf('NRMSE: \t\t\t%f\nvMag Error: \t\t%f\nAbsolute Angle Error: \t%f\n\n',vNRMSE_RBF,vMagErr_RBF,angErr_RBF);
pause(1);

%% Plot all
% Plotting results from all methods.

figure,imshow_flowmag(...
    cat(2,imMask,imMask,imMask,imMask,imMask,imMask),...
    cat(2,vx,vxN,vxDFWsm,vxDFWsms,vxFDM,vxRBF),...
    cat(2,vy,vyN,vyDFWsm,vyDFWsms,vyFDM,vyRBF),...
    cat(2,vz,vzN,vzDFWsm,vzDFWsms,vzFDM,vzRBF),...
    vMax,[1,2,3]);

title('Original,Noisy,DFW,DFW/spin,FDM and RBF','FontSize',14)
pause(1)

%% Visualize DivFree Wavelet Coefficients

% Wavelet coefficients for original
minSize = 8*ones(1,3);
[vx,vy,vz] = maskIM(imMask,vx,vy,vz);
[wcdf1,wcdf2,wcn,numLevels,wcSizes] = dfwavelet_forward(vx,vy,vz,minSize,res);
wcdf1_tile = wcTile(wcdf1,numLevels,wcSizes);
wcdf2_tile = wcTile(wcdf2,numLevels,wcSizes);
wcn_tile = wcTile(wcn,numLevels,wcSizes);

% Wavelet coefficients for noisy
[wcdf1N,wcdf2N,wcnN,numLevels,wcSizes] = dfwavelet_forward(vxN,vyN,vzN,minSize,res);
wcdf1N_tile = wcTile(wcdf1N,numLevels,wcSizes);
wcdf2N_tile = wcTile(wcdf2N,numLevels,wcSizes);
wcnN_tile = wcTile(wcnN,numLevels,wcSizes);


% Wavelet coefficeints for div free wavelet w/ cycle spinning
[wcdf1sms,wcdf2sms,wcnsms,numLevels,wcSizes] = dfwavelet_forward(vxDFWsms,vyDFWsms,vzDFWsms,minSize,res);
wcdf1sms_tile = wcTile(wcdf1sms,numLevels,wcSizes);
wcdf2sms_tile = wcTile(wcdf2sms,numLevels,wcSizes);
wcnsms_tile = wcTile(wcnsms,numLevels,wcSizes);   


% Plot
figure,imshow3f(abs(cat(2,wcn_tile,wcnN_tile,wcnsms_tile)),[1,2,3],round(N/4)),title('Non-divfree wavelet subbands for Original, Noisy and DFW','FontSize',14);




