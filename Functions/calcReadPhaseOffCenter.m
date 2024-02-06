function [phase_off_c read_off_c] = calcReadPhaseOffCenter(prot)
% Aaron Hess
% University of Oxford
% Calculate phase and read off center
% prot is structure generated from ascconvtrans

    TRANSVERSE = 1; SAGITTAL = 2; CORONAL = 3;
    [mfilepath,NAME,EXT] = fileparts(mfilename('fullpath')) ;
    addpath([mfilepath, filesep, 'ascii_prot']);

%     prot = ascconvtrans(fname);
    
    pos(1) = str2double(prot.sSliceArray.asSlice(1).sPosition.dSag);
    pos(2) = str2double(prot.sSliceArray.asSlice(1).sPosition.dCor);
    pos(3) = str2double(prot.sSliceArray.asSlice(1).sPosition.dTra);

    norm(1) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dSag);
    norm(2) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dCor);
    norm(3) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dTra);
    inPlaneRoation = str2double(prot.sSliceArray.asSlice(1).dInPlaneRot);
        
    [phase read orienation] = fGSLCalcPRS(norm, inPlaneRoation);
    % Slice pos.PEvec = phase offcenter
    phase_off_c = pos*phase';
    read_off_c = pos*read';
    
    rmpath([mfilepath, filesep, 'ascii_prot']);
end