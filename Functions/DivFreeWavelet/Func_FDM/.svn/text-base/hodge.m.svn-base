function [vxDF,vyDF,vzDF,vxCF,vyCF,vzCF] = hodge(vx,vy,vz,res)
%
% [vxDF,vyDF,vzDF,vxCF,vyCF,vzCF] = hodge(vx,vy,vz,res)
%
% Hodge decomposition using the fast Poisson solver
%
% Inputs:
%     vx, vy, vz
%     res     -   resolution
%
% Outputs:
%     vxDF,vyDF,vzDF    -  DivFree
%     vxCF,vyCF,vzCF    -  Non-divfree
%
% (c) Frank Ong 2013

if (isempty(res))
    res = [1,1,1];
end

% d = div(v) = div(grad(u) = laplacian(u)
d = div(vx,vy,vz,res);
% u = laplacian^-1(u)
u = poisson(d,res);
% vCF = grad(u)
[vxCF,vyCF,vzCF] = grad(u,res);

% vDF = v - vCF
vxDF = vx-vxCF;
vyDF = vy-vyCF;
vzDF = vz-vzCF;
