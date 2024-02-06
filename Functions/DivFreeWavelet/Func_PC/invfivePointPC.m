function [vx,vy,vz,imMag] = invfivePointPC(s1,s2,s3,s4,s5,M1)
%
% [vx,vy,vz,imMag] = invfivePointPC(s1,s2,s3,s4,s5,M1)
%
% The following implements 5-point phase contrast encoding as described in:
% Kevin M. Johnson and Michael Markl
% Improved SNR in Phase Contrast Velocimetry With Five-Point Balanced Flow Encoding
% 
% Inputs: 
%     s1,s2,s3,s4,s5 - velocity encodes
%     M1      -   First moment
%
% Outputs:
%     vx,vy,vz-   flow data
%     imMag   -   image magnitude
%
%       
% (c) Frank Ong 2013

gamma = 26753;

A1 = [-1,-1,-1;
    1,1,-1;
    1,-1,1;
    -1,1,1;];
A0 = [A1;
    0,0,0];
Ani = A1*pinv(A0);
A1inv = pinv(A1);

ph1 = angle(s1);
ph2 = angle(s2);
ph3 = angle(s3);
ph4 = angle(s4);
ph5 = angle(s5);

ni1 = 2*pi*round(1/(2*pi)*(Ani(1,1)*ph1+Ani(1,2)*ph2+Ani(1,3)*ph3+Ani(1,4)*ph4+Ani(1,5)*ph5-ph1));
ni2 = 2*pi*round(1/(2*pi)*(Ani(2,1)*ph1+Ani(2,2)*ph2+Ani(2,3)*ph3+Ani(2,4)*ph4+Ani(2,5)*ph5-ph2));
ni3 = 2*pi*round(1/(2*pi)*(Ani(3,1)*ph1+Ani(3,2)*ph2+Ani(3,3)*ph3+Ani(3,4)*ph4+Ani(3,5)*ph5-ph3));
ni4 = 2*pi*round(1/(2*pi)*(Ani(4,1)*ph1+Ani(4,2)*ph2+Ani(4,3)*ph3+Ani(4,4)*ph4+Ani(4,5)*ph5-ph4));

vx = 1/(gamma*M1)*(A1inv(1,1)*(ph1+ni1)+A1inv(1,2)*(ph2+ni2)+A1inv(1,3)*(ph3+ni3)+A1inv(1,4)*(ph4+ni4));
vy = 1/(gamma*M1)*(A1inv(2,1)*(ph1+ni1)+A1inv(2,2)*(ph2+ni2)+A1inv(2,3)*(ph3+ni3)+A1inv(2,4)*(ph4+ni4));
vz = 1/(gamma*M1)*(A1inv(3,1)*(ph1+ni1)+A1inv(3,2)*(ph2+ni2)+A1inv(3,3)*(ph3+ni3)+A1inv(3,4)*(ph4+ni4));
imMag = (abs(s1)+abs(s2)+abs(s3)+abs(s4)+abs(s5))/5;