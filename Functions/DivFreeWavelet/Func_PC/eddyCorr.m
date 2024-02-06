function [vxEddy,vyEddy,vzEddy] = eddyCorr(vx,vy,vz,calib)
%
% [vxEddy,vyEddy,vzEddy] = eddyCorr(vx,vy,vz,calib)
%
% 2nd order polynomial fitting for eddy current correction
% 
% Inputs: 
%     vx,vy,vz -   Input velocity data
%     calib    -   Calibration region (static flow)
%
% Outputs:
%     vxEddy,vyEddy,vzEddy - Correction terms
%
%
% Example:
%     [vxEddy,vyEddy,vzEddy] = eddyCorr(vx,vy,vz,calib)
%     vx = vx - vxEddy;
%     vy = vy - vyEddy;
%     vz = vz - vzEddy;
%       
% (c) Frank Ong 2013
%


[xGrid,yGrid,zGrid] = ndgrid(calib(1,1):calib(1,2),calib(2,1):calib(2,2),calib(3,1):calib(3,2));
x = xGrid(:);
y = yGrid(:);
z = zGrid(:);
A = [ones(numel(xGrid),1),x,y,z,x.^2,y.^2,z.^2,x.*y,x.*z,y.*z];
[xBig,yBig,zBig] = ndgrid(1:size(vx,1),1:size(vx,2),1:size(vx,3));

%% Eddy-current fitting

vxCalib = vx(calib(1,1):calib(1,2),calib(2,1):calib(2,2),calib(3,1):calib(3,2));
poly = A\(vxCalib(:));
vxEddy = poly(1)+poly(2)*xBig+poly(3)*yBig+poly(4)*zBig+...
poly(5)*xBig.^2+poly(6)*yBig.^2+poly(7)*zBig.^2+...
poly(8)*xBig.*yBig+poly(9)*xBig.*zBig+poly(10)*yBig.*zBig;

%%
vyCalib = vy(calib(1,1):calib(1,2),calib(2,1):calib(2,2),calib(3,1):calib(3,2));
poly = A\(vyCalib(:));
vyEddy = poly(1)+poly(2)*xBig+poly(3)*yBig+poly(4)*zBig+...
poly(5)*xBig.^2+poly(6)*yBig.^2+poly(7)*zBig.^2+...
poly(8)*xBig.*yBig+poly(9)*xBig.*zBig+poly(10)*yBig.*zBig;

%%
vzCalib = vz(calib(1,1):calib(1,2),calib(2,1):calib(2,2),calib(3,1):calib(3,2));
poly = A\(vzCalib(:));
vzEddy = poly(1)+poly(2)*xBig+poly(3)*yBig+poly(4)*zBig+...
poly(5)*xBig.^2+poly(6)*yBig.^2+poly(7)*zBig.^2+...
poly(8)*xBig.*yBig+poly(9)*xBig.*zBig+poly(10)*yBig.*zBig;
