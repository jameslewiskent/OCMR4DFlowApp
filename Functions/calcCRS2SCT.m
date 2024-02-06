function CRS2SCT = calcCRS2SCT(info)
% Aaron Hess
% Modified from Coded by Alper Yaman, Feb 2009
% info = dicominfo(fname);

    ipp = info.ImagePositionPatient;
    iop = info.ImageOrientationPatient;
    ps  = info.PixelSpacing;

    Tipp=[1 0 0 ipp(1); 0 1 0 ipp(2); 0 0 1 ipp(3); 0 0 0 1];

    row = iop(1:3);  
    col = iop(4:6); 
    slc = cross(row,col);

    %R = [r(1) c(1) s(1) 0; r(2) c(2) s(2) 0; r(3) c(3) s(3) 0; 0 0 0 1];
    R = diag([1 1 1 1]);
    R(1:3,1:3) = [row,col,slc];

    if info.MRAcquisitionType=='3D' % 3D turboflash
        S = [ps(1) 0 0 0; 0 ps(2) 0 0; 0 0 info.SliceThickness 0 ; 0 0 0 1];
    else % 2D epi dti
        if(isfield(info,'SpacingBetweenSlices'))
            S = [ps(1) 0 0 0;0 ps(2) 0 0;0 0 info.SpacingBetweenSlices 0;0 0 0 1];
        else
            S = [ps(1) 0 0 0; 0 ps(2) 0 0; 0 0 info.SliceThickness 0 ; 0 0 0 1];
        end
    end
    T0 = [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    
    T0(1,1) = 0;
    T0(1,2) = 1;
    T0(2,1) = 1;
    T0(2,2) = 0;

    CRS2SCT = Tipp * R * S * T0;
    
end