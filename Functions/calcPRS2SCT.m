function [prs2sct orienation inPlaneRoation] = calcPRS2SCT(prot, fType)
% Aaron HEss
% University of Oxford
% August 2012
% Extract orientation of Phase read slice from ascconv section of dicom
% fType = 1; for imaging (defalt)
% fType = 2; for spectroscopy
% orientation is a Siemens deffinition
% 
% prot is from ascconvtrans(fname);

    TRANSVERSE = 1; SAGITTAL = 2; CORONAL = 3;

    if(~exist('fType','var'))
        fType = 1;
    end
    
%     prot = ascconvtrans(fname);

    if (fType==2)
        matrix_phase = str2double(prot.sSpecPara.lFinalMatrixSizePhase);
        matrix_read  = str2double(prot.sSpecPara.lFinalMatrixSizeRead);
        matrix_par = 1;
        phase_res = 1;
    else
        matrix_phase = str2double(prot.sKSpace.lPhaseEncodingLines);
        matrix_read = str2double(prot.sKSpace.lBaseResolution);
        matrix_par = str2double(prot.sKSpace.lImagesPerSlab);
        phase_res = str2double(prot.sKSpace.dPhaseResolution);

    end

    if (fType<2)
        %pos(1) = str2double(prot.sSliceArray.asSlice(1).sPosition.dSag);
        %% ATH 06/2013 PRS to SCT should have the same reference location
        %pos(2) = str2double(prot.sSliceArray.asSlice(1).sPosition.dCor);
        %pos(3) = str2double(prot.sSliceArray.asSlice(1).sPosition.dTra);
        pos = [0 0 0];
        
        dim(1) = str2double(prot.sSliceArray.asSlice(1).dPhaseFOV)/(matrix_phase/phase_res);
        dim(2) = str2double(prot.sSliceArray.asSlice(1).dReadoutFOV)/matrix_read;
        dim(3) = str2double(prot.sSliceArray.asSlice(1).dThickness)/matrix_par;

        norm(1) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dSag);
        norm(2) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dCor);
        norm(3) =  str2double(prot.sSliceArray.asSlice(1).sNormal.dTra);
        inPlaneRoation = str2double(prot.sSliceArray.asSlice(1).dInPlaneRot);

    else
        pos(1) = str2double(prot.sSpecPara.sVoI.sPosition.dSag);
        pos(2) = str2double(prot.sSpecPara.sVoI.sPosition.dCor);
        pos(3) = str2double(prot.sSpecPara.sVoI.sPosition.dTra);
        % NB DIFFERENT FOR CSI
        if(fType == 3)
            dim(1) = str2double(prot.sSpecPara.sVoI.dPhaseFOV)/(matrix_phase/phase_res);
            dim(2) = str2double(prot.sSpecPara.sVoI.dReadoutFOV)/matrix_read;
            dim(3) = str2double(prot.sSpecPara.sVoI.dThickness)/matrix_par;
        else
            dim(1) = str2double(prot.sSliceArray.asSlice(1).dPhaseFOV)/(matrix_phase/phase_res);
            dim(2) = str2double(prot.sSliceArray.asSlice(1).dReadoutFOV)/matrix_read;
            dim(3) = str2double(prot.sSliceArray.asSlice(1).dThickness)/matrix_par;
        end
        norm(1) =  str2double(prot.sSpecPara.sVoI.sNormal.dSag);
        norm(2) =  str2double(prot.sSpecPara.sVoI.sNormal.dCor);
        norm(3) =  str2double(prot.sSpecPara.sVoI.sNormal.dTra);
        inPlaneRoation = str2double(prot.sSpecPara.sVoI.dInPlaneRot);

    end

    [phase read orienation] = fGSLCalcPRS(norm, inPlaneRoation);
    prs2sct = zeros(4,4);
    prs2sct(1:3,1:3) = [phase(:), read(:), norm(:)];
    prs2sct(:,4) = [pos(:); 1];
  
end

