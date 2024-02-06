function [s1,s2,s3,s4,s5] = fivePointPC(vx,vy,vz,imMag,ph0,M1)
%
% [s1,s2,s3,s4,s5] = fivePointPC(vx,vy,vz,imMag,ph0,M1)
%
% The following implements 5-point phase contrast encoding as described in:
% Kevin M. Johnson and Michael Markl
% Improved SNR in Phase Contrast Velocimetry With Five-Point Balanced Flow Encoding
% 
% Inputs: 
%     vx,vy,vz-   flow data
%     imMag   -   image magnitude
%     ph0     -   reference phase
%     M1      -   First moment
%
% Outputs:
%     s1,s2,s3,s4,s5 - velocity encodes
%       
% (c) Frank Ong 2013
gamma = 26753;

A1 = [-1,-1,-1;
    1,1,-1;
    1,-1,1;
    -1,1,1;];
A0 = [A1;
    0,0,0];


s1 = imMag.*exp(1i*((A0(1,1)*vx+A0(1,2)*vy+A0(1,3)*vz)*gamma*M1+ph0));
s2 = imMag.*exp(1i*((A0(2,1)*vx+A0(2,2)*vy+A0(2,3)*vz)*gamma*M1+ph0));
s3 = imMag.*exp(1i*((A0(3,1)*vx+A0(3,2)*vy+A0(3,3)*vz)*gamma*M1+ph0));
s4 = imMag.*exp(1i*((A0(4,1)*vx+A0(4,2)*vy+A0(4,3)*vz)*gamma*M1+ph0));
s5 = imMag.*exp(1i*((A0(5,1)*vx+A0(5,2)*vy+A0(5,3)*vz)*gamma*M1+ph0));
