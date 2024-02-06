function MovieWriter(movieName,img,scale)
    % Render out movie frames
    sz = size(img);
    if(length(sz)>3)
        NPHS = sz(4);
        % mosaic third dimenssion
        img_show = Mozaic(img);
    else
        img_show = abs(img);
        NPHS = sz(3);
    end
    
    if(~exist('scale','var'))
        scale = [0 prctile(abs(img(:)),98)];
    end
    
    img_show = 255*(img_show-scale(1))/(scale(2)-scale(1));
    
    
    
    clear mv; frameCounter = 1; % Re-initalize
    fg = figure;
    subplot('position', [0 0 1 1]) 
    set(gcf, 'color', 'black');
    for ii = 1:NPHS
        imshow(imresize(img_show(:,:,ii),4),[0 255],'border','tight');
        axis image; axis off; colormap gray; 
        
        %mv(frameCounter) = getframe(gcf);
        mv(frameCounter).cdata = repmat(uint8(imresize(img_show(:,:,ii),4/2.15)),[1 1 3]);
        mv(frameCounter).colormap =[];
        
        frameCounter = frameCounter + 1;
    end

    
    if(~isempty(movieName))
        % Save to AVI file
        writerObj = VideoWriter(movieName,'MPEG-4');
        writerObj.FrameRate = NPHS;
        writerObj.Quality = 100;
        % writerObj.profile = 'MPEG-4';
        open(writerObj);
        writeVideo(writerObj,mv);
        close(writerObj);
    end
    close(fg)
end