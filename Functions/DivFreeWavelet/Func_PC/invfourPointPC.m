function [vx,vy,vz,ph0,imMag] = invfourPointPC(s1,s2,s3,s4,M1)
%
% [vx,vy,vz,ph0,imMag] = invfourPointPC(s1,s2,s3,s4,M1)
%
% The following implements 4-point phase contrast encoding
%
% Inputs: 
%     s1,s2,s3,s4 - velocity encodes
%     M1      -   First moment
%
% Outputs:
%     vx,vy,vz-   flow data
%     imMag   -   image magnitude
%       
% (c) Frank Ong 2013

gamma = 26753;


% A = [-1    -1    -1     +1;
%     +1    +1    -1     +1;
%     +1    -1    +1     +1;
%     -1    +1    +1     +1];


% A = [-1    -1    -1     1;
%     1    -1    -1     1;
%     1     1    -1     1;
%     1     1     1     1];


vx = 1/(gamma*M1)*angle(conj(s1).*s2)/2;
vy = 1/(gamma*M1)*angle(conj(s2).*s3)/2;
vz = 1/(gamma*M1)*angle(conj(s3).*s4)/2;
ph0 = angle(s1.*s4)/2;
imMag = (abs(s1)+abs(s2)+abs(s3)+abs(s4))/4;