% Save figures or output images for solved displacement and strain fields

%%
% Find img name
[~,imgname,imgext] = fileparts([fileNameLeft{2,ImgSeqNum},'\',fileNameLeft{1,ImgSeqNum}]);
%%
if isempty(DICpara.outputFilePath)
    DICpara.outputFilePath = uigetdir;
    outputVariables = {'DispU','DispV','DispW','exx','exy','eyy','eyy2'};
    for i = 1:length(outputVariables)
        tempfolder_path = [DICpara.outputFilePath,'\',outputVariables{i}];
        if ~exist(tempfolder_path, 'dir')  
            mkdir(tempfolder_path);  
        end
    end
end

if DICpara.MethodToSaveFig == 1
    %% jpg
    figure(1); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    print([DICpara.outputFilePath,'\DispU\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'],'-djpeg','-r300');

    figure(2); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    print([DICpara.outputFilePath,'\DispV\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'],'-djpeg','-r300');

    figure(3); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim([-0.025,0.025]); end
    print([DICpara.outputFilePath,'\DispW\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispW'],'-djpeg','-r300')

    % % figure(4); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim([-0.025,0.025]); end
    % % print([DICpara.outputFilePath,'\DispA\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispA'],'-djpeg','-r300')

    figure(4); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim([-0.015,0.015]); end
    print([DICpara.outputFilePath,'\exx\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx'],'-djpeg','-r300')

    figure(5); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([0,0.025]); end
    print([DICpara.outputFilePath,'\exy\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy'],'-djpeg','-r300')

    figure(6); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([-0.025,0]); end
    print([DICpara.outputFilePath,'\eyy\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy'],'-djpeg','-r300')

    figure(7); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([-0.025,0]); end
    print([DICpara.outputFilePath,'\eyy2\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy_2'],'-djpeg','-r300')

    % 
    % % figure(8); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([0,0.07]); end
    % % print([DICpara.outputFilePath,'\DispU\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxprincipal'],'-djpeg','-r300')
    % % 
    % % figure(9); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([0,0.07]); end
    % % print([DICpara.outputFilePath,'\DispU\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises'],'-djpeg','-r300')
    % % 
    % figure(7); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([0,0.07]); end
    % print([DICpara.outputFilePath,'\dwdx\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_dwdx'],'-djpeg','-r300')
    % 
    % figure(8); if DICpara.OrigDICImgTransparency == 0, colormap jet;  clim([0,0.07]); end
    % print([DICpara.outputFilePath,'\dwdy\',imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_dwdy'],'-djpeg','-r300')

elseif DICpara.MethodToSaveFig == 2
    %% pdf
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'];
    figure(3); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'];
    figure(4); if DICpara.OrigDICImgTransparency == 0, colormap jet; clim auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx'];
    figure(5); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy'];
    figure(6); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy'];
    figure(7); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([-0.015,0.015]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_max'];
    figure(8); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([0,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_min'];
    figure(9); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([-0.025,0]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxshear'];
    figure(10); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises'];
    figure(11); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); clim([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    
else
    %% fig
    fprintf('Please modify codes manually in Section 8.');
    figure(3); colormap(coolwarm(128)); clim auto; 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU.fig']);
    figure(4); colormap(coolwarm(128)); clim auto; 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV.fig']);
    figure(5); colormap(coolwarm(128)); clim([-0.05,0.1]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx.fig']);
    figure(6); colormap(coolwarm(128)); clim([-0.05,0.05]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy.fig']);
    figure(7); colormap(coolwarm(128)); clim([-0.1,0.05]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy.fig']);
    figure(8); colormap(coolwarm(128)); clim([0,0.025]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_max.fig']);
    figure(9); colormap(coolwarm(128)); clim([-0.025,0]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_min.fig']);
    figure(10); colormap(coolwarm(128)); clim([0,0.07]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxshear.fig']);
    figure(11); colormap(coolwarm(128)); clim([0,0.07]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises.fig']);
    
    
end





