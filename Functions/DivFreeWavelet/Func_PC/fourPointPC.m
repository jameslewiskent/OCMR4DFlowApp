function [s1,s2,s3,s4] = fourPointPC(vx,vy,vz,imMag,ph0,M1)
%
% [s1,s2,s3,s4] = fourPointPC(vx,vy,vz,imMag,ph0,M1)
%
% The following implements 4-point phase contrast encoding
%
% Inputs: 
%     vx,vy,vz-   flow data
%     imMag   -   image magnitude
%     ph0     -   reference phase
%     M1      -   First moment
%
% Outputs:
%     s1,s2,s3,s4 - velocity encodes
%       
% (c) Frank Ong 2013

gamma = 26753;


% A = [-1    -1    -1     1;
%     1    -1    -1     1;
%     1     1    -1     1;
%     1     1     1     1];


s1 = imMag.*exp(1i*((-vx-vy-vz)*gamma*M1+ph0));
s2 = imMag.*exp(1i*((vx-vy-vz)*gamma*M1+ph0));
s3 = imMag.*exp(1i*((vx+vy-vz)*gamma*M1+ph0));
s4 = imMag.*exp(1i*((vx+vy+vz)*gamma*M1+ph0));
