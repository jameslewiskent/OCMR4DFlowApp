function app = GetListOfContours(app, user_input)
            % Find list of contours in selected directory
            contour_file_extension = '.cvi42wsx';
            list = dir(fullfile(app.directoryPath,['*',contour_file_extension]));
            
            % No countours found, give error
            if isempty(list)
                uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'Could not find any contour files in specified directory.','An error occured');
                app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Could not find any contour files in specified directory. '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
                error('Could not find any contour files in specified directory.')
            end
            
            % Update contour drop down            
            for contour_n = 1:size(list,1)
            contour_list{1,contour_n} = list(contour_n).name;
            end
            
            % If more than one, and not a group directory add option 'All
            % of above'
            if ~app.GroupCheckBox.Value && size(list,1) > 1 && user_input == 1
                contour_list{1,size(list,1) + 1} = '-Analyse all of above-';
            end
            
            app.ContourDropDown.Items = contour_list;
            app.ContourDropDown.ItemsData = contour_list;
            app.ContourDropDown.Value = contour_list(1);
end

