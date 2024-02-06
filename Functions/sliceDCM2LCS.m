function img_out = sliceDCM2LCS(img_in, orientation, isPESwapped)
% Aaron Hess
% University of Oxford
% Reorient slice object from Dicom orientation to logical (phase read slc)
% Rotate by 90 / -270 is Shift 0,0 to m,0
% Rotate by 270 / -90 is Shift 0,0 to 0,n
% Updated May 2014 for multi dimensional input

    TRANSVERSE = 1; SAGITTAL = 2; CORONAL = 3;
    
    img_out = [];
    [d1 d2 d3 d4 d5] = size(img_in);
    
    
    for i3 = 1:d3
        for i4 = 1:d4
            for i5 = 1:d5
                img = img_in(:,:,i3,i4,i5);
                switch(orientation)

                    case SAGITTAL
                        if(~isPESwapped)
                            % rot -90 or +270
                            img_LCS = rot90(img,3);  %rot by 270
                        else
                            % do nothing
                            img_LCS = img;
                        end

                    case CORONAL
                        if(~isPESwapped)
                            % mirr COL
                            % rot -270 or +90
                            img_LCS = mirror_COL(img);
                            img_LCS = rot90(img_LCS);
                        else
                            % mirr LIN
                            img_LCS = mirror_LIN(img);
                        end

                    case TRANSVERSE
                        if(~isPESwapped)
                            %mirr COL
                            %img_LCS = img;
                            img_LCS = mirror_COL(img);
                        else
                            % mirr COL
                            % rot -270 or +90
                            img_LCS = mirror_COL(img);
                            img_LCS = rot90(img_LCS);
                        end
                    otherwise
                        error('sliceDCM2LCS:: unknown orientation')
                end

                img_LCS = mirror_LIN(img_LCS);  % Do not know why i need to do this...
                
                if isempty(img_out)
                    img_out = zeros(size(img_LCS,1),size(img_LCS,2), d3, d4,d5);
                end
                img_out(:,:,i3,i4,i5) = img_LCS;
            end
        end
    end
                
end

function img_m = mirror_COL(img)
    img_m = img(:,end:-1:1);
end

function img_m = mirror_LIN(img)
    img_m = img(end:-1:1,:);
end