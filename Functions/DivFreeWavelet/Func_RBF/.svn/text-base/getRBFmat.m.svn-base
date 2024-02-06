function rbfM = getRBFmat(x,y,z,dist,res)
%
% rbfKer = getRBFmat(x,y,z,dist,res)
%
% Generates RBF matrix for each spacial point
% dist is a parameter that controls the size of the kernel
% 
%
% Inputs:
%     x,y,z   -   point
%     dist    -   sort of the distance to midpoint, controls the support
%     res     -   Resolution [1x3]
%
% Outputs:
%     rbfKer  -   structure that contains the kernels
%
% (c) Frank Ong 2013

    r = [x*res(1);y*res(2);z*res(3)];
    rNorm = norm(r);
    rbfM = ((1-rNorm^2/(2*dist^2))*eye(3)+1/(2*dist^2)*(r*r'))*exp(-rNorm^2/(2*dist^2));
end