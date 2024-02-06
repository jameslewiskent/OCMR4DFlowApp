function [vNRMSE,vMagErr,angErr] = calcVelError(imMask,vx_ref,vy_ref,vz_ref,vx_com,vy_com,vz_com)
%
% [vNRMSE,vMagErr,angErr] = calcVelError(imMask,vx_ref,vy_ref,vz_ref,vx_com,vy_com,vz_com)
%
% Calculates errors in three different metrics
% 
% Inputs: 
%     imMask               -   mask that defines region of interest
%     vx_ref,vy_ref,vz_ref -   reference velocities
%     vx_com,vy_com,vz_com -   compared velocities
%
% Outputs:
%     vNRMSE               -   normalized root mean squared error
%     vMagErr              -   normalized root mean squared error for
%                              velocity magnitude
%     angErr               -   direction error
%       
% (c) Frank Ong 2013

imMask = imMask>0;
vx_ref = vx_ref(imMask);
vy_ref = vy_ref(imMask);
vz_ref = vz_ref(imMask);
vx_com = vx_com(imMask);
vy_com = vy_com(imMask);
vz_com = vz_com(imMask);

% Normalized Root Mean Square Err for velocity magnitude (speed)
vMag_ref = sqrt(vx_ref.^2+vy_ref.^2+vz_ref.^2);
vMag_com = sqrt(vx_com.^2+vy_com.^2+vz_com.^2);
vMax = max(vMag_ref(:));
vMagErr = sqrt(mean((vMag_ref-vMag_com).^2))/vMax;

% Normalized Root Mean Square Err
vNRMSE = sqrt(mean((vx_ref-vx_com).^2+(vy_ref-vy_com).^2+(vz_ref-vz_com).^2))/vMax;

% Angle Error
angErr = abs(vx_ref.*vx_com + vy_ref.*vy_com + vz_ref.*vz_com );
denom = vMag_ref.*vMag_com;
angErr = 1-angErr./denom;
angErr(vMag_ref==0) = 0;
angErr(vMag_com==0) = 0;
angErr = mean(angErr);