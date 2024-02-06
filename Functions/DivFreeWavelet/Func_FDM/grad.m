function [fx,fy,fz] = grad(in,res)
% [fx,fy,fz] = grad(in,res)
% calculates gradient (adjoint to div)

if (isempty(res))
    res = [1,1,1];
end
fx = (circshift(in,[-1,0,0])-in)/(res(1));
fy = (circshift(in,[0,-1,0])-in)/(res(2));
fz = (circshift(in,[0,0,-1])-in)/(res(3));
 