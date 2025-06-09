% This script is used for % Set
folderPath = 'C:\Users\Zixiang Tong\OneDrive - The University of Texas at Austin\Documents\MATLABCodes\2D_ALDIC-master\CalibrationImages\CheckBoard\SimCheckerBoardCal\16-mm\R'; % 替换为你的文件夹路径

% 获取文件夹中所有的.tif文件
tifFiles = dir(fullfile(folderPath, '*.tiff'));

% 循环处理每个.tif文件，将其重命名为.tiff
for i = 1:length(tifFiles)
    oldFileName = fullfile(folderPath, tifFiles(i).name); % 获取旧文件名（含路径）
    newFileName = fullfile(folderPath, [tifFiles(i).name(1:end-5), '.tif']); % 构建新文件名（含路径）
    movefile(oldFileName, newFileName); % 重命名文件
end

disp('批量重命名完成。');
