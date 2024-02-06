% Create a protocol structure with some manditory null data
% the asc files do not generate an entry if a variable is null so therefore
% initialization is nesseary
function prot = prot_empty()
    prot.sAdjData.sAdjVolume.sPosition.dSag = {'0'};
    prot.sAdjData.sAdjVolume.sPosition.dCor = {'0'};
    prot.sAdjData.sAdjVolume.sPosition.dTra = {'0'};

    prot.sAdjData.sAdjVolume.sNormal.dSag = {'0'};
    prot.sAdjData.sAdjVolume.sNormal.dCor = {'0'};
    prot.sAdjData.sAdjVolume.sNormal.dTra = {'0'};
    prot.sAdjData.sAdjVolume.dInPlaneRot = {'0'};
    
    prot.sSliceArray.asSlice(1).sPosition.dSag = {'0'};
    prot.sSliceArray.asSlice(1).sPosition.dCor = {'0'};
    prot.sSliceArray.asSlice(1).sPosition.dTra = {'0'};
    
    prot.sSliceArray.asSlice(1).sNormal.dSag = {'0'};
    prot.sSliceArray.asSlice(1).sNormal.dCor = {'0'};
    prot.sSliceArray.asSlice(1).sNormal.dTra = {'0'};
    
    prot.sSliceArray.asSlice(1).dInPlaneRot =  {'0'};
end