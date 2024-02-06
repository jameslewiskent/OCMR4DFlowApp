function out = convfft3(in,ker)
%
% out = convfft3(in,ker)
%
% Convolution using FFT. Fast when kernel is big
% 
% Inputs: 
%     in      -   Input Data
%     ker     -   Kernel
%
% Outputs:
%     out     -  output
%       
% (c) Frank Ong 2013

newSize = size(in)+size(ker)+1;

out2 = ifftn( fftn(in,newSize) .* fftn(ker,newSize) );
outInd = ceil((size(ker)-1)/2);
out = out2(outInd(1) + (1:size(in,1)) , outInd(2) + (1:size(in,2)) , outInd(3) + (1:size(in,3)));