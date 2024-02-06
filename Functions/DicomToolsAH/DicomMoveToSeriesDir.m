function DicomMoveToSeriesDir(dir_in,app)
% Aaron Hess
% University of Oxford 2014
% Move dicoms into folders named by their series name
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Arranging Dicom folder structure');

[dicomData] = processDicomDir(dir_in, '*'); % Do not do recursive!
fclose('all');
for Study_n = 1:size(dicomData.study,2) % Number of studies
    for iSer = 1:length(dicomData.study(Study_n).series)
        series = dicomData.study(Study_n).series(iSer);
        new_dir_name = [dir_in filesep num2str(series.SeriesNumber), '_', series.SeriesDescription];
        [isGood, msg, ~] = mkdir(new_dir_name);
        if(~isGood)
            app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + [' - Could not create directory ',new_dir_name, ' because ',msg]); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        else
            for ifile = 1:length(series.instance)
                [~, fname, fext] = fileparts(series.instance(ifile).Filename);
                % JK for reasons I do not(!) understand, *sometimes* this movefile fails
                % and it tries to move a file FROM the destination?!
                % error message 'cannot find file specified'
                [isGood, msg, ~] = movefile(series.instance(ifile).Filename,[new_dir_name filesep fname fext],'f'); % JK added 'f'
                if(~isGood)
                    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + [' - Could not move file ',fname, ' because ',msg]); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
                    
                    try % Try doing a system level move file if the above hasn't worked
                        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Attempting a system level move file.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
                        if    ispc
                            system(['copy ' series.instance(ifile).Filename ' ' [new_dir_name filesep fname fext] ])
                        else
                            system(['cp -r ' series.instance(ifile).Filename ' ' [new_dir_name filesep fname fext] ])
                        end
                    catch
                        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - System level move file also failed. Check dicoms are the expected format.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
                    end
                end
                progbar.Value = (((iSer-1)*length(series.instance))+ifile)/(length(dicomData.study(Study_n).series)*length(series.instance));
            end
        end
    end
end
progbar.Value = 1;
close(progbar);
end