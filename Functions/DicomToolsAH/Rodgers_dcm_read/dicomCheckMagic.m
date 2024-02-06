function [bIsDicom] = dicomCheckMagic(strFilename)
    thisDir = dir(strFilename);
    
    if numel(thisDir) ~= 1 || thisDir.isdir || thisDir.bytes < 132 % 128 + 4
        bIsDicom = false;
        return
    end
    
    fid = fopen(strFilename,'r');
    if fid == -1
        warning('ShMOLLI:dicomCheckMagic','Cannot open file "%s".',strFilename)
        bIsDicom = false;
        return
    end
    c = onCleanup(@() fclose(fid));
    
    fseek(fid, 128, 'bof');
    
    magic = fread(fid,4,'*char*1');
    
    bIsDicom = strcmp(reshape(magic,1,[]),'DICM');
end
