function [vxRBF,vyRBF,vzRBF] = rbfDenoise(vx,vy,vz,imMask,radius,res,nIter,doplot)
%
% [vxRBF,vyRBF,vzRBF] = rbfDenoise(vx,vy,vz,imMask,radius,res,nIter,doplot)
%
% The following implements divergence-free radial basis function denoising 
% as described in:
% Busch J, Kozerke S., et al. MRM 2012
% Construction of divergence-free velocity fields from cine 3D phase-contrast flow measurements. 
%
%
% Inputs:
%     vx,vy,vz  -   3d matrices of velocities
%     imMask    -   uncertainty function
%     radius    -   radius of kernel, kernel has size of (2*radius+1)^3
%     res       -   resolution
%     nIter     -   number of iteration of LSQR
%     doplot    -   plot during iteration. Slower when enabled.
%
% Outputs:
%     vxRBF,vyRBF,vzRBF
%
% (c) Frank Ong 2013

fprintf('RBF Denoising...\n');

FOV = size(vx);

vx = vx.*imMask;
vy = vy.*imMask;
vz = vz.*imMask;

vMax = sqrt(max(vx(:).^2+vy(:).^2+vz(:).^2));

v = [vx(:);vy(:);vz(:)];

L = length(v);

rbfKer = getRBFker(radius,res);

coeff = zeros(L,1);

if doplot
    % doplot does gradient descent instead of lsqr
    for i = 1:nIter
        
        [coeff,flag] = lsqr( @(x,transp_flag) RBFconv(x,transp_flag,imMask,FOV,rbfKer),...
            v,[],1,[],[],coeff);
        
        vi = RBFconv(coeff,'notransp',imMask,FOV,rbfKer);
        
        vxRBF = reshape(vi(1:L/3),FOV);
        vyRBF = reshape(vi((L/3+1):(2*L/3)),FOV);
        vzRBF = reshape(vi((2*L/3+1):end),FOV);
        
        figure(76),clf,imshow_flowmag(imMask,vxRBF,vyRBF,vzRBF,vMax,[1,2,3],round(FOV(3)/2))
        title(sprintf('RBF Iter: %d',i),'FontSize',14), drawnow
        
    end
    close(76)
else
    [coeff,flag] = lsqr( @(x,transp_flag) RBFconv(x,transp_flag,imMask,FOV,rbfKer),...
        v,[],nIter,[],[],coeff);
end


v = RBFconv(coeff,'notransp',imMask,FOV,rbfKer);

vxRBF = reshape(v(1:L/3),FOV);
vyRBF = reshape(v((L/3+1):(2*L/3)),FOV);
vzRBF = reshape(v((2*L/3+1):end),FOV);



fprintf('done\n');