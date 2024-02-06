function u=poisson(f,res)
%
% u = poisson ( f, res)
%
% First order Poisson solver using the FFT with periodic extension
% Solves \laplacian u = f
%
% Inputs:
%     f       -   div of flow field
%     res     -   Resolution [1x3]
%
% Outputs:
%     u       -   flux of flow field
%
% (c) Frank Ong 2013


if (isempty(res))
    res = [1,1,1];
end

[nx,ny,nz]=size(f);
[X,Y,Z] = ndgrid(0:nx-1,0:ny-1,0:nz-1);

% lambdas are Fourier Transform of Discrete Laplacian: [-2,1,0...0,1];
lambdax= -4*(sin(X*pi/nx)).^2/res(1)^2;
lambday= -4*(sin(Y*pi/ny)).^2/res(2)^2;
lambdaz= -4*(sin(Z*pi/nz)).^2/res(3)^2;
mu = ( lambdax + lambday + lambdaz);
mu(1) = 1;

u = fftn(f)./mu;
u(1) = 0;
u = real(ifftn(u));