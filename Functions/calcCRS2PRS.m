function [crs2prs, orienation, isPESwaped] = calcCRS2PRS(prot, info)
% Aaron Hess
% Calculate transformation matrix from column row slice to phase read slice
% University of oxford
% August 2012
% prot is from ascconvtrans(fname)
% info is from dicominfo

    crs2sct = calcCRS2SCT(info);

    [prs2sct, orienation, inPlaneRoation] = calcPRS2SCT(prot);

    crs2prs = prs2sct\crs2sct;

    % Calculate if PE is swapped
    if(inPlaneRoation > pi/4) && (inPlaneRoation < 3*pi/4)
        isPESwaped = true;
    else
        isPESwaped = false;
    end

% % Determin if any axes need to be reflected
% TRANSVERSE = 1; SAGITTAL = 2; CORONAL = 3;
% switch(orienation)
%     case TRANSVERSE
%     case SAGITTAL
%         % reflect phase encode
%         crs2prs(1,:) = -crs2prs(1,:);
%     case CORONAL
%         if(isPESwaped)
%             crs2prs(2,:) = -crs2prs(2,:);
%         else
%             crs2prs(2,:) = -crs2prs(1,:);
%         end
% end
end