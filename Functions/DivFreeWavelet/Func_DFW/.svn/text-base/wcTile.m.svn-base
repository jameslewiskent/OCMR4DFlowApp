function out = wcTile(in,numLevels,wcSizes)
%
% out = wcTile(in,numLevels,wcSizes)
%
% Reshape 3D wavelet coefficients into matrix and normalize each subband
% 
% Inputs: 
%     in          -   input wavelet coefficients
%     numLevels   -   Number of wavelet decomposition levels
%     wc_sizes    -   wavelet subband sizes 
%                     (The first three numbers correspond to the scaling
%                     subband, the second group corresponds to the smallest
%                     detail subband and it gets bigger and bigger
%
% Outputs:
%     out         -   tiled wavelet coefficients
%
%
%       
% (c) Frank Ong 2013

outdims = wcSizes(1:3)*2^numLevels;

out = zeros(outdims);

ind = 1;
for l = 1:numLevels
    blockDim = wcSizes((3*l+1):(3*l+3));
    blockSize = prod(blockDim);
    
    x0 = 1:blockDim(1);
    y0 = 1:blockDim(2);
    z0 = 1:blockDim(3);
    x1 = wcSizes(1)*2^(l-1) + (1:blockDim(1));
    y1 = wcSizes(2)*2^(l-1) + (1:blockDim(2));
    z1 = wcSizes(3)*2^(l-1) + (1:blockDim(3));
    
    if (l==1)
        out(x0,y0,z0) = reshape(in(ind:(ind+blockSize-1)),blockDim);
        out(x0,y0,z0) = normalize(out(x0,y0,z0));
        ind = ind+blockSize;
    end
    
    out(x1,y0,z0) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x1,y0,z0) = normalize(out(x1,y0,z0));
    ind = ind+blockSize;
    
    out(x0,y1,z0) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x0,y1,z0) = normalize(out(x0,y1,z0));
    ind = ind+blockSize;
    
    out(x1,y1,z0) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x1,y1,z0) = normalize(out(x1,y1,z0));
    ind = ind+blockSize;
    
    out(x0,y0,z1) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x0,y0,z1) = normalize(out(x0,y0,z1));
    ind = ind+blockSize;
    
    out(x1,y0,z1) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x1,y0,z1) = normalize(out(x1,y0,z1));
    ind = ind+blockSize;
    
    out(x0,y1,z1) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x0,y1,z1) = normalize(out(x0,y1,z1));
    ind = ind+blockSize;
    
    out(x1,y1,z1) = reshape(in(ind:(ind+blockSize-1)),blockDim);
    out(x1,y1,z1) = normalize(out(x1,y1,z1));
    ind = ind+blockSize;
end

function out = normalize(in)
m = max(abs(in(:)));
if m~=0
    out = in/m;
else
    out = 0;
end