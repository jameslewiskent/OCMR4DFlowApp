function [vxRec,vyRec,vzRec] = fdmDenoise(vx,vy,vz,res)
%
% [vxRec,vyRec,vzRec] = fdmDenoise(vx,vy,vz,res)
%
% The following implements finite difference method denoising 
% as described in:
% Song SM, Pelc NJ., et al. JMRI 1993
% Noise reduction in three-dimensional phase-contrast MR velocity measurements.
%
%
% Inputs:
%     vx,vy,vz  -   3d matrices of velocities
%     res       -   resolution
%
% Outputs:
%     vxRec,vyRec,vzRec
%
% (c) Frank Ong 2013

if (isempty(res))
    res = [1,1,1];
end
[vxRec,vyRec,vzRec,a,b,c] = hodge(vx,vy,vz,res);