function y = RBFconv(x,transp_flag,imMask,FOV,rbfKer)
%
% y = RBFconv(x,transp_flag,imMask,FOV,rbfKer)
%
% RBF convolution function for LSQR
%
% Inputs:
%     x         -   input data ( [vx;vy;vz] or [coef1;coef2;coef3] )
%     transp_flag-  transpose flag. Since RBF is symmetric, the only
%                   difference is in applying imMask
%     imMask    -   uncertainty function
%     FOV       -   image size
%     rbfKer    -   RBF kernels
%
% Outputs:
%     t         -   output data ( [coef1;coef2;coef3] or [vx;vy;vz]  )
%
% (c) Frank Ong 2013

L = length(x);

x1 = reshape(x(1:L/3),FOV);
x2 = reshape(x((L/3+1):(2*L/3)),FOV);
x3 = reshape(x((2*L/3+1):end),FOV);

if strcmp(transp_flag,'transp')
    % Velocity to RBF coefficients
    x1 = x1.*imMask;
    x2 = x2.*imMask;
    x3 = x3.*imMask;
end

if (size(rbfKer.rbf11,1) < 9)
    y11 = convn(x1,rbfKer.rbf11,'same');
    y12 = convn(x2,rbfKer.rbf12,'same');
    y13 = convn(x3,rbfKer.rbf13,'same');
    
    y21 = convn(x1,rbfKer.rbf21,'same');
    y22 = convn(x2,rbfKer.rbf22,'same');
    y23 = convn(x3,rbfKer.rbf23,'same');
    
    y31 = convn(x1,rbfKer.rbf31,'same');
    y32 = convn(x2,rbfKer.rbf32,'same');
    y33 = convn(x3,rbfKer.rbf33,'same');
else
    
    y11 = convfft3(x1,rbfKer.rbf11);
    y12 = convfft3(x2,rbfKer.rbf12);
    y13 = convfft3(x3,rbfKer.rbf13);
    
    y21 = convfft3(x1,rbfKer.rbf21);
    y22 = convfft3(x2,rbfKer.rbf22);
    y23 = convfft3(x3,rbfKer.rbf23);
    
    y31 = convfft3(x1,rbfKer.rbf31);
    y32 = convfft3(x2,rbfKer.rbf32);
    y33 = convfft3(x3,rbfKer.rbf33);
end

y1 = y11 + y12 + y13;
y2 = y21 + y22 + y23;
y3 = y31 + y32 + y33;


if strcmp(transp_flag,'notransp')
    y1 = y1.*imMask;
    y2 = y2.*imMask;
    y3 = y3.*imMask;
end


y = [y1(:);y2(:);y3(:)];
