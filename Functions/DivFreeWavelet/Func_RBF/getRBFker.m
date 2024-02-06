function rbfKer = getRBFker(radius,res)
%
% rbfKer = getRBFker(radius,res)
%
% Generates RBF convolution kernel for given radius
% The size of the kernel is (2*radius+1)^3
%
% Inputs:
%     r       -   radius of kernel
%     res     -   Resolution [1x3]
%
% Outputs:
%     rbfKer  -   structure that contains the kernels
%
% (c) Frank Ong 2013

rbf11 = zeros(2*radius+1,2*radius+1,2*radius+1);
rbf12 = rbf11;
rbf13 = rbf11;

rbf21 = rbf11;
rbf22 = rbf11;
rbf23 = rbf11;

rbf31 = rbf11;
rbf32 = rbf11;
rbf33 = rbf11;

for z = (-radius:radius)
    for y = (-radius:radius)
        for x = (-radius:radius)
            i = 1+radius+x;
            j = 1+radius+y;
            k = 1+radius+z;
            
            rbfM = getRBFmat(x,y,z,min(radius*res)/2,res);

            % rbfM is actually symmetric but for code readability, we save 
            % all of the kernels
            rbf11(i,j,k) = rbfM(1,1);
            rbf12(i,j,k) = rbfM(1,2);
            rbf13(i,j,k) = rbfM(1,3);

            rbf21(i,j,k) = rbfM(2,1);
            rbf22(i,j,k) = rbfM(2,2);
            rbf23(i,j,k) = rbfM(2,3);
                
            rbf31(i,j,k) = rbfM(3,1);
            rbf32(i,j,k) = rbfM(3,2);
            rbf33(i,j,k) = rbfM(3,3);
                
        end
    end
end

rbfKer.rbf11 = rbf11;
rbfKer.rbf12 = rbf12;
rbfKer.rbf13 = rbf13;
rbfKer.rbf21 = rbf21;
rbfKer.rbf22 = rbf22;
rbfKer.rbf23 = rbf23;
rbfKer.rbf31 = rbf31;
rbfKer.rbf32 = rbf32;
rbfKer.rbf33 = rbf33;


