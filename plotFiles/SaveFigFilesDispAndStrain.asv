function SaveFigFilesDispAndStrain(DICpara, fileNameLeft, ImgSeqNum, varargin)
%SAVEFIGFILESDISPANDSTRAIN Save displacement and strain figures to various file formats
%
% USAGE:
%   SaveFigFilesDispAndStrain(DICpara, fileNameLeft, ImgSeqNum)
%   SaveFigFilesDispAndStrain(DICpara, fileNameLeft, ImgSeqNum, 'FileFormat', 'jpg')
%   SaveFigFilesDispAndStrain(DICpara, fileNameLeft, ImgSeqNum, 'Resolution', 600)
%   SaveFigFilesDispAndStrain(DICpara, fileNameLeft, ImgSeqNum, 'FigureList', [1,2,3])
%
% INPUTS:
%   DICpara     - DIC parameters structure
%   fileNameLeft- Left camera file names cell array
%   ImgSeqNum   - Current image sequence number
%
% OPTIONAL PARAMETERS:
%   'FileFormat'  - Output file format: 'jpg'(default), 'png', 'pdf', 'eps', 'fig', 'tiff'
%   'Resolution'  - Output resolution in DPI (default: 300)
%   'FigureList'  - List of figure numbers to save (default: all open figures)
%   'OutputPath'  - Custom output path (default: uses DICpara.outputFilePath or prompts user)
%   'Prefix'      - Custom filename prefix (default: uses image name)
%   'CreateSubfolders' - Create subfolders for each variable (default: true)
%
% AUTHOR: Modified for 3D-Stereo-ALDIC by Zach Tong

%% Parse input arguments
p = inputParser;
addRequired(p, 'DICpara', @isstruct);
addRequired(p, 'fileNameLeft', @iscell);
addRequired(p, 'ImgSeqNum', @isnumeric);
addParameter(p, 'FileFormat', 'jpg', @(x) ismember(lower(x), {'jpg', 'jpeg', 'png', 'pdf', 'eps', 'fig', 'tiff', 'tif'}));
addParameter(p, 'Resolution', 300, @isnumeric);
addParameter(p, 'FigureList', [], @isnumeric);
addParameter(p, 'OutputPath', '', @ischar);
addParameter(p, 'Prefix', '', @ischar);
addParameter(p, 'CreateSubfolders', true, @islogical);

parse(p, DICpara, fileNameLeft, ImgSeqNum, varargin{:});

fileFormat = lower(p.Results.FileFormat);
resolution = p.Results.Resolution;
figureList = p.Results.FigureList;
outputPath = p.Results.OutputPath;
prefix = p.Results.Prefix;
createSubfolders = p.Results.CreateSubfolders;

%% Generate image name and prefix
try
    [~, imgname, ~] = fileparts([fileNameLeft{2,ImgSeqNum}, '\', fileNameLeft{1,ImgSeqNum}]);
catch
    try
        [~, imgname, ~] = fileparts(fileNameLeft{1,ImgSeqNum});
    catch
        imgname = sprintf('Image_%04d', ImgSeqNum);
    end
end

if isempty(prefix)
    prefix = imgname;
end

%% Set up output path and folders
if isempty(outputPath)
    if isfield(DICpara, 'outputFilePath') && ~isempty(DICpara.outputFilePath)
        outputPath = DICpara.outputFilePath;
    else
        outputPath = uigetdir(pwd, 'Select output directory for figures');
        if outputPath == 0
            fprintf('Save operation cancelled by user.\n');
            return;
        end
        DICpara.outputFilePath = outputPath;
    end
end

% Define output variables and their corresponding folder names
outputVariables = {'DispU', 'DispV', 'DispW', 'DispMag', ...
                   'exx', 'exy', 'eyy', 'eyz', 'exz', 'ezz', ...
                   'dwdx', 'dwdy', 'principal_max', 'principal_min', ...
                   'max_shear', 'vonMises', 'AllFigures'};

% Create subfolders if requested
if createSubfolders
    for i = 1:length(outputVariables)
        tempfolder_path = fullfile(outputPath, outputVariables{i});
        if ~exist(tempfolder_path, 'dir')
            mkdir(tempfolder_path);
        end
    end
end

%% Get list of figures to save
if isempty(figureList)
    figureList = findobj('Type', 'figure');
    figureList = [figureList.Number];
end

if isempty(figureList)
    fprintf('No figures found to save.\n');
    return;
end

%% Set format-specific parameters
switch fileFormat
    case {'jpg', 'jpeg'}
        formatFlag = '-djpeg';
        fileExt = '.jpg';
    case 'png'
        formatFlag = '-dpng';
        fileExt = '.png';
    case 'pdf'
        formatFlag = '-dpdf';
        fileExt = '.pdf';
    case 'eps'
        formatFlag = '-depsc';
        fileExt = '.eps';
    case {'tiff', 'tif'}
        formatFlag = '-dtiff';
        fileExt = '.tif';
    case 'fig'
        formatFlag = '';
        fileExt = '.fig';
    otherwise
        error('Unsupported file format: %s', fileFormat);
end

%% Save figures
fprintf('Saving figures in %s format...\n', upper(fileFormat));

% Define figure-to-variable mapping
figureMapping = containers.Map('KeyType', 'int32', 'ValueType', 'char');
figureMapping(1) = 'DispU';
figureMapping(2) = 'DispV';
figureMapping(3) = 'DispW';
figureMapping(4) = 'DispMag';
figureMapping(5) = 'exx';
figureMapping(6) = 'exy';
figureMapping(7) = 'eyy';
figureMapping(8) = 'eyz';
figureMapping(9) = 'exz';
figureMapping(10) = 'ezz';
figureMapping(11) = 'dwdx';
figureMapping(12) = 'dwdy';
figureMapping(13) = 'principal_max';
figureMapping(14) = 'principal_min';
figureMapping(15) = 'max_shear';
figureMapping(16) = 'vonMises';

for figNum = figureList
    try
        figure(figNum);
        
        % Generate filename components
        if isfield(DICpara, 'winsize') && isfield(DICpara, 'winstepsize')
            paramStr = sprintf('_WS%d_ST%d', DICpara.winsize, DICpara.winstepsize);
        else
            paramStr = '';
        end
        
        % Determine variable name and subfolder
        if figureMapping.isKey(figNum)
            varName = figureMapping(figNum);
        else
            varName = sprintf('Figure_%d', figNum);
        end
        
        % Generate full filename
        if createSubfolders
            if ismember(varName, outputVariables)
                fullPath = fullfile(outputPath, varName, [prefix, paramStr, '_', varName, fileExt]);
            else
                fullPath = fullfile(outputPath, 'AllFigures', [prefix, paramStr, '_', varName, fileExt]);
            end
        else
            fullPath = fullfile(outputPath, [prefix, paramStr, '_', varName, fileExt]);
        end
        
        % Apply colormap and axis settings if needed
        if isfield(DICpara, 'OrigDICImgTransparency') && DICpara.OrigDICImgTransparency == 0
            colormap('turbo'); % Use modern colormap
            clim('auto');
        end
        
        % Save the figure
        if strcmp(fileFormat, 'fig')
            savefig(fullPath);
        else
            resolutionStr = sprintf('-r%d', resolution);
            print(fullPath, formatFlag, resolutionStr);
        end
        
        fprintf('Saved: %s\n', fullPath);
        
    catch ME
        fprintf('Error saving figure %d: %s\n', figNum, ME.message);
    end
end

fprintf('Figure saving completed. Files saved to: %s\n', outputPath);

end





