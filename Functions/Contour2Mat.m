function [out_data, out_name] = Contour2Mat(d_path,contour_file, infoNative, info_out,aux_name)
% Aaron Hess
% University of Oxford
% April 2018
% Read in a contour file from CVI42 
% convert to scanner coordinates?
% create a mask according to a specific data set
% Extract a set of contours, form a mask for largest and smallest VOI
% if aux_name is given the contour file is amended with that name

%% 10. Read in contours and generate masks
% match each image to a contour
% NOTE NOT ASSUMING ALL CONTOURS ARE NEEDED
% Must strip all contours not matching our label...
% iscontours = (~isempty(contour_files))&&iscine;
% if(~( iscontours))
%     fprintf('No cine or contours provided\n');
% else

    if(~exist('aux_name','var'))
        aux_name = '';
    end
    %dk only
    SliceL = uint8(infoNative.Private_0051_100d);

    [~, sl_index] = UpdateList(infoNative.order.sliceList,SliceL);
    CRS2SCT_native = calcCRS2SCT(infoNative);  % Note something not quite right here...
    CRS2SCT_native(:,4) = CRS2SCT_native*([-1 -1 -sl_index 1])'; % correct for slice which info belongs to, and so that index[1 1 1] is what was [0 0 0]
    CRS2SCT_native = CRS2SCT_native/1000; % convert to m (from mm)
    CRS2SCT_native(4,4) = 1;
    
    % these are contours based on our destination geometry
    if(~exist('info_out','var'))%||~exist('CRS2SCT_out','var'))
        info_out = infoNative;
        CRS2SCT_out = CRS2SCT_native;
    else
        [~, sl_index] = UpdateList(info_out.order.sliceList,info_out.SliceL);
        CRS2SCT_out = calcCRS2SCT(info_out);  % Note something not quite right here...
        CRS2SCT_out(:,4) = CRS2SCT_out*([-1 -1 -sl_index 1])'; % correct for slice which info belongs to, and so that index[1 1 1] is what was [0 0 0]
        CRS2SCT_out = CRS2SCT_out/1000; % convert to m (from mm)
        CRS2SCT_out(4,4) = 1;        
    end

    mat_file_cont = fullfile(d_path,[extractBefore(contour_file,'.'),'.mat']);
 

    fprintf('loading contours\n');

    labels = {};
    no_uid = cellfun(@isempty,infoNative.uids);
    infoNative.uids(no_uid) = {'empty'};

    [contours, ~] = ExtractCMR42Contours(fullfile(d_path,contour_file));

        % add label to this to allow multiple labels
        for iImg = 1:length(contours)
            for iCont = 1:length(contours(iImg).str)
                [labels,~] = UpdateList(labels,contours(iImg).str(iCont).Label);
            end
        end
        if length(labels)>1
            [iLabel,~] = listdlg('ListString',labels,'PromptString','select which contour label to use', 'SelectionMode','single');
        else
            iLabel  = 1;
        end
        [~,c_nameBase,~] = fileparts(contour_file);

%         for iLabel = 1:length(labels)
            c_name = [c_nameBase,labels{iLabel}];
        %     cont_match = zeros(length(contours),8);
            cont_vol = zeros(size(infoNative.uids));
            areaToVol = 1e-3*prod(diag(CRS2SCT_native))/4/4;  % units ml (4 as ,Points are x 4 ...)

        %     phsCine2Flow = NTime/NPhsCine;

            BWNative = zeros(infoNative.size);

            clear('pt_all');
            pt_all{infoNative.size(6)} = [];
            dcm_found = false;

            for iImg = 1:length(contours)
                is_match = strfind(infoNative.uids,contours(iImg).uid);
                ind = find(~cellfun(@isempty,is_match));
                % create a label list too...
                if(isempty(ind))
                    warning(['Could not match contours to cine, contour file: ',contour_file]);
                    continue;
                end
                [a, b, c, d, e, f, g] = ind2sub(size(infoNative.uids),ind);  % a is slice d is time frame
        %         cont_match(iImg,2:8) = [a, b, c, d, e, f,g];  
        %         cont_match(iImg,1) = ind;
                for iCont = 1:length(contours(iImg).str)
                    if(isempty(contours(iImg).str(iCont).Label))
                        contours(iImg).str(iCont).Label = 'na';
                    end                    
                    if(strcmp(labels{iLabel},contours(iImg).str(iCont).Label))
                        cont_vol(a,b,c,d,e,f,g) = areaToVol*polyarea(contours(iImg).str(iCont).Points(1,:),contours(iImg).str(iCont).Points(2,:));


                        % convert each contour into 4D flow coordinates then create a mask
                        pt = [contours(iImg).str(iCont).Points(2:-1:1,:)/4; repmat(a,[1,size(contours(iImg).str(iCont).Points,2)]); repmat(1,[1,size(contours(iImg).str(iCont).Points,2)])];
                        %ptFlow = cine2flow*pt;
                        BWNative(:,:,a,b,c,d,e,f,g) = roipoly(BWNative(:,:,1),pt(2,:),pt(1,:));
                        % now match ECG time stamp?
                        pt_all{d} = [pt_all{d},pt];
                        dcm_found = true;
                    end
                end
            end

%             if(~dcm_found)
%                 warning(['Image match not found for contours labeled: ',labels{iLabel}]);
%                 continue;
%             end
            % remove contrast, te,ind, dimenssion
            if(size(BWNative,8)>1)
                BWNative = max(BWNative,[],8);
            end
            if(size(BWNative,4)>1)
                BWNative = max(BWNative,[],4);
            end
            if(size(BWNative,5)>1)
                BWNative = max(BWNative,[],5);
            end
            if(size(BWNative,7)>1)
                BWNative = max(BWNative,[],7);
            end
            if(size(BWNative,9)>1)
                BWNative = max(BWNative,[],9);
            end
            

            BWNative = squeeze(BWNative);
            
            volBW = sum(sum(sum(BWNative,1),2),3);

            
            % remove time slices with no data
            is_time = find(volBW>0);
            [~,i]=min(volBW(is_time));
            i_min = is_time(i);
            [~,i]=max(volBW(is_time));
            i_max = is_time(i);
            
            if(length(is_time)<2)
                %fprintf('\nlength(is_time) = %g\n',length(is_time));
                warning('Contour2Mat Error: less than two contours found to match data')
            end
            
            
            % create a flow mask from contours
            BWMax = zeros(info_out.size(1:3));
            BWMin = zeros(info_out.size(1:3));
            
            [X, Y, Z, T, C] = ind2sub(size(BWMin),1:numel(BWMin));

            % Tcine = T/phsCine2Flow;  %this streaches cine contours to fit flow
%             Tcine = T;
%             Tcine(T>size(BWNative,4))=size(BWNative,4); % this assumes same timing just shorter cycle
            inCine = CRS2SCT_native\CRS2SCT_out*[X;Y;Z;C];
            inCine(1,inCine(1,:)<1) = 1;
            inCine(2,inCine(2,:)<1) = 1;
            inCine(3,inCine(3,:)<1) = 1;
            inCine(1,inCine(1,:)>size(BWNative,1)) = size(BWNative,1);
            inCine(2,inCine(2,:)>size(BWNative,2)) = size(BWNative,2);
            inCine(3,inCine(3,:)>size(BWNative,3)) = size(BWNative,3);

            inCine = round(inCine);
            

            BWNativeMin = BWNative(:,:,:,i_min);
            BWNativeMax = BWNative(:,:,:,i_max);
            BWMax(1:end) = BWNativeMax(sub2ind(infoNative.size(1:3),inCine(1,:),inCine(2,:),inCine(3,:)));
            BWMin(1:end) = BWNativeMin(sub2ind(infoNative.size(1:3),inCine(1,:),inCine(2,:),inCine(3,:)));
% 
%             cont.BWCine = BWNative;
%             cont.BWFlow = BWOut;
%             cont.pt_all = pt_all;
% 
%             cont.vol_t = vol_t;


%             var.([c_name,'Mask']) = repmat(BWFlow,[1 1 1 rep_cycle]);


    save(mat_file_cont,'BWMax','BWMin','BWNative','CRS2SCT_native','pt_all');
    
    out_data.BWMax=BWMax;
    out_data.BWMin=BWMin;
    out_data.BWNative = BWNative;
    out_data.pt_all = pt_all;
    out_data.CRS2SCT_native = CRS2SCT_native;
    
    out_name = mat_file_cont;
end  % funtion
