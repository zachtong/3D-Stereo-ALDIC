% Save figures or output images for solved displacement and strain fields

%%
% Find img name
[~,imgname,imgext] = fileparts([fileNameLeft{2,ImgSeqNum},'\',fileNameLeft{1,ImgSeqNum}]);
%%
if isempty(DICpara.outputFilePath)
    DICpara.outputFilePath = uigetdir;
    outputVariables = {'DispU','DispV','DispW'};
    for i = 1:length(outputVariables)
        tempfolder_path = [DICpara.outputFilePath,'\',outputVariables{i}];
        if ~exist(tempfolder_path, 'dir')  
            mkdir(tempfolder_path);  
        end
    end
end

    %% jpg
    figure(1); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    print([DICpara.outputFilePath,'\DispU\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'],'-djpeg','-r300');

    figure(2); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    print([DICpara.outputFilePath,'\DispV\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'],'-djpeg','-r300');

    figure(3); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim([-0.025,0.025]); end
    print([DICpara.outputFilePath,'\DispW\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispW'],'-djpeg','-r300')

  
 




