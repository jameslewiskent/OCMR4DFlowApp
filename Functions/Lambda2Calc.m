function lambda2 = Lambda2Calc(flow,spacing)
% Aaron Hess
% University of Oxford
% November 2019
% Calculate Lambda2 for vortex identification in a 4D flow flow field
% flow is arranged as Nx x Ny x Nz x Nt x 3
% spacing is voxel spacing in x y and z

    X = 1; Y=2; Z=3;
    [Nx, Ny, Nz, Nt, v]  = size(flow);

    % clculate gradient,
    % NOTE THIS IS CRUDE, it should be  a sensile interpolation
    dxyz = zeros(v,v,Nx,Ny,Nz,Nt);

    for it = 1:Nt
        for id = 1:3
            [tx,ty,tz] = imgradientxyz(flow(:,:,:,it,id));
            dxyz(id,X,:,:,:,it) = tx/spacing(1);
            dxyz(id,Y,:,:,:,it) = ty/spacing(2);
            dxyz(id,Z,:,:,:,it) = tz/spacing(3);
        end
    end


    [S, ~,  SdS, OdO] = calculateSymmetry(reshape(dxyz,[9,Nx*Ny*Nz*Nt]));

    eig3d = arrayfun(@(ind) eig(SdS(:,:,ind)+OdO(:,:,ind)), 1:size(S,3), 'uniformOutput', false);
    eigv = cat(3, eig3d{:});
    lambda2 = reshape(eigv(2,1,:),[Nx,Ny,Nz,Nt]);

end


% Function Adapted from TurbTools.m 
% http://turbulence.pha.jhu.edu/matlabanalysistools.aspx
 % Function to calculate symmetric- and antisymmetric parts of
% tensor, and their squares. Input is a 9 by N matrix, where the 9
% rows represent the 9 components of the tensor, and N the number
% of points for which we have this tensor
function [S, O,  SdS, OdO] =  calculateSymmetry(gradient)
    Jt = reshape(gradient, 3, 3, length(gradient));
    J = permute(Jt, [2 1 3]);
    S = (J+Jt)/2;
    O = (J-Jt)/2;
    St = permute(S, [2 1 3]);
    Ot = permute(O, [2 1 3]);

%     multiply3dS = arrayfun(@(ind) S(:, :, ind) * St(:, :, ind), 1:size(S,3), 'uniformOutput', false);
%     multiply3dO = arrayfun(@(ind) O(:, :, ind) * Ot(:, :, ind), 1:size(S,3), 'uniformOutput', false);
%     SSt = cat(3, multiply3dS{:});
%     OOt = cat(3, multiply3dO{:});

    SdS = zeros(3,3,length(S));
    OdO = zeros(3,3,length(O));
    for i = 1:3
        for k = 1:3
            for j = 1:3
                SdS(i,k,:) = SdS(i,k,:) + S(i,j,:) .* S(j,k,:);
                OdO(i,k,:) = OdO(i,k,:) + O(i,j,:) .* O(j,k,:);
            end
        end
    end
end