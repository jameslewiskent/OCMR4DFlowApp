function [vx,vy,vz,imMag] = genCouette(vMag,N,Rout,Rin,Win)
% [vx,vy,vz,imMag] = genCouette(vMag,N,Rout,Rin,Win)
%
% Generates a Couette flow data with given parameters
%
% see http://en.wikipedia.org/wiki/Couette_flow

FOV = [N,N,N];
n = Rin/Rout;
A = -Win*n^2/(1-n^2);
B = Win*Rin^2/(1-n^2);

vx = zeros(FOV);
vy = zeros(FOV);
vz = zeros(FOV);
imMag = zeros(FOV);

center = (N+1)/2;
scale = max(abs(A*Rout+B/Rout),abs(A*Rin+B/Rin));

for k = (1:N)
    for j = (1:N)
        for i = (1:N)
            x = i-center;
            y = j-center;
            r = sqrt(x^2+y^2);
            if (r<Rout) && (r>Rin)
                theta = atan2(y,x);
                vx(i,j,k) = -vMag*(A*r+B/r)*sin(theta)/scale;
                vy(i,j,k) = vMag*(A*r+B/r)*cos(theta)/scale;
                imMag(i,j,k) = 1;
            end
        end
    end
end