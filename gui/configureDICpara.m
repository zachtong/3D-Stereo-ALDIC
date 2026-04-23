function DICpara = configureDICpara(DICparaIn)
%CONFIGUREDICPARA  Interactive GUI for filling a DICpara struct.
%
%   DICpara = configureDICpara()              % start from defaults
%   DICpara = configureDICpara(DICparaIn)     % pre-fill from struct
%
%   Opens a modal uifigure with four sections: Data / DIC parameters /
%   Strain & smoothing / Output. User picks folders, enters numeric
%   values, selects options, and clicks Run. On Run, validates via
%   validateDICpara and returns the filled struct. On Cancel, returns
%   [] (caller should bail out).

if nargin < 1, DICparaIn = setDICparaDefaults(); end

% State container passed to all callbacks
state.DICpara = DICparaIn;
state.accepted = false;

% ---------- uifigure ----------
fig = uifigure('Name', '3D-Stereo-ALDIC Configuration', ...
               'Position', [200, 100, 720, 780], ...
               'Resize', 'off', ...
               'CloseRequestFcn', @(~,~) onCancel());

mainGrid = uigridlayout(fig, [6, 1]);
mainGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 40};
mainGrid.RowSpacing = 6;
mainGrid.Padding = [10 10 10 10];

% ---------- 1. Data panel ----------
[~, dataHandles] = buildDataPanel(mainGrid, state.DICpara);

% ---------- 2. DIC parameters panel ----------
[~, dicHandles] = buildDICPanel(mainGrid, state.DICpara);

% ---------- 3. Strain & smoothing panel ----------
[~, strainHandles] = buildStrainPanel(mainGrid, state.DICpara);

% ---------- 4. Output panel ----------
[~, outputHandles] = buildOutputPanel(mainGrid, state.DICpara);

% ---------- Status line ----------
statusLabel = uilabel(mainGrid, 'Text', '', 'FontColor', [0.6 0 0], ...
                     'HorizontalAlignment', 'left', 'WordWrap', 'on');

% ---------- Button row ----------
btnGrid = uigridlayout(mainGrid, [1, 3]);
btnGrid.ColumnWidth = {'1x', 100, 100};
cancelBtn = uibutton(btnGrid, 'push', 'Text', 'Cancel', ...
                    'ButtonPushedFcn', @(~,~) onCancel());
cancelBtn.Layout.Column = 2;
runBtn = uibutton(btnGrid, 'push', 'Text', 'Run ▶', ...
                 'BackgroundColor', [0.2 0.6 0.2], ...
                 'FontColor', 'white', 'FontWeight', 'bold', ...
                 'ButtonPushedFcn', @(~,~) onRun());
runBtn.Layout.Column = 3;

% ---------- Wait modal ----------
uiwait(fig);

if state.accepted
    DICpara = state.DICpara;
else
    DICpara = [];
end


% ================== Callbacks (nested functions) ==================

    function onCancel()
        state.accepted = false;
        if isvalid(fig), delete(fig); end
    end

    function onRun()
        try
            collectValues();
            validateDICpara(state.DICpara);
            validatePaths();
            state.accepted = true;
            if isvalid(fig), delete(fig); end
        catch err
            statusLabel.Text = ['ERROR: ' err.message];
            drawnow;
        end
    end

    function collectValues()
        % Data
        state.DICpara.imgFolderLeft   = dataHandles.imgLeftEdit.Value;
        state.DICpara.imgFolderRight  = dataHandles.imgRightEdit.Value;
        state.DICpara.maskFolder      = dataHandles.maskEdit.Value;
        state.DICpara.calibrationFile = dataHandles.calibEdit.Value;
        state.DICpara.calibrationMethod = dataHandles.calibFormatDD.Value;
        state.DICpara.DICIncOrNot     = dataHandles.modeDD.Value;

        % DIC
        state.DICpara.winsize            = dicHandles.winsize.Value;
        state.DICpara.winstepsize        = dicHandles.winstepsize.Value;
        state.DICpara.winsizeMin         = dicHandles.winsizeMin.Value;
        state.DICpara.NewFFTSearchDistance = [dicHandles.fftDistX.Value, dicHandles.fftDistY.Value];
        state.DICpara.ClusterNo          = dicHandles.clusterNo.Value;
        state.DICpara.UseGlobal          = dicHandles.useGlobal.Value;

        % Strain
        state.DICpara.strain_size    = strainHandles.strainSize.Value;
        state.DICpara.Smooth2DTimes  = strainHandles.smooth2D.Value;
        state.DICpara.Smooth3DTimes  = strainHandles.smooth3D.Value;
        state.DICpara.transformDisp  = double(strainHandles.transformDisp.Value);

        % Output
        state.DICpara.PlotMode         = outputHandles.plotMode.Value;
        state.DICpara.Image2PlotResults = outputHandles.image2Plot.Value;
        state.DICpara.MethodToSaveFig  = outputHandles.figFormat.Value;
        state.DICpara.outputFilePath   = outputHandles.outputPath.Value;
    end

    function validatePaths()
        p = state.DICpara;
        if isempty(p.imgFolderLeft)  || ~exist(p.imgFolderLeft,  'dir'), error('Left images folder not set or missing.'); end
        if isempty(p.imgFolderRight) || ~exist(p.imgFolderRight, 'dir'), error('Right images folder not set or missing.'); end
        if isempty(p.maskFolder)     || ~exist(p.maskFolder,     'dir'), error('Mask folder not set or missing.'); end
        if isempty(p.calibrationFile)|| ~exist(p.calibrationFile,'file'),error('Calibration file not set or missing.'); end
    end
end


% ================== Panel builders ==================

function [panel, h] = buildDataPanel(parent, p)
panel = uipanel(parent, 'Title', '1. Data', 'FontWeight', 'bold');
g = uigridlayout(panel, [6, 3]);
g.ColumnWidth = {200, '1x', 80};
g.RowHeight = repmat({28}, 1, 6);
g.RowSpacing = 4;

addRow(g, 1, 'Left images folder:', @() pickFolder('Select LEFT images folder'));
h.imgLeftEdit = g.Children(end-1);
h.imgLeftEdit.Value = safeStr(p, 'imgFolderLeft');
h.imgLeftEdit.Tooltip = 'Folder containing left-camera image frames (common formats).';

addRow(g, 2, 'Right images folder:', @() pickFolder('Select RIGHT images folder'));
h.imgRightEdit = g.Children(end-1);
h.imgRightEdit.Value = safeStr(p, 'imgFolderRight');
h.imgRightEdit.Tooltip = 'Folder containing right-camera image frames.';

addRow(g, 3, 'Mask folder:', @() pickFolder('Select mask folder'));
h.maskEdit = g.Children(end-1);
h.maskEdit.Value = safeStr(p, 'maskFolder');
h.maskEdit.Tooltip = 'Masks named *_0.* (left) and *_1.* (right). Acc mode uses the first pair only; inc mode needs one pair per frame.';

addRow(g, 4, 'Calibration file:', @() pickFile('Select calibration file'));
h.calibEdit = g.Children(end-1);
h.calibEdit.Value = safeStr(p, 'calibrationFile');
h.calibEdit.Tooltip = 'Stereo calibration file. Format must match the ''Calibration format'' dropdown.';

uilabel(g, 'Text', 'Calibration format:');
h.calibFormatDD = uidropdown(g, ...
    'Items', {'MATLAB Calibrator (case 0)', 'MatchID (case 1)', 'MCC (case 2)', 'DICe XML (case 3)', 'OpenCorr (case 4)'}, ...
    'ItemsData', [0 1 2 3 4], ...
    'Value', fallbackVal(p, 'calibrationMethod', 3));
h.calibFormatDD.Tooltip = 'Choose the calibration file format.';
uilabel(g, 'Text', '');
h.calibFormatDD.Layout.Row = 5;
h.calibFormatDD.Layout.Column = 2;

uilabel(g, 'Text', 'DIC mode:');
h.modeDD = uidropdown(g, ...
    'Items', {'Accumulative (small disp)', 'Incremental (large disp)'}, ...
    'ItemsData', [0 1], ...
    'Value', fallbackVal(p, 'DICIncOrNot', 0));
h.modeDD.Tooltip = 'Accumulative: first frame is the reference for all. Incremental: each frame referenced to the previous.';
uilabel(g, 'Text', '');
h.modeDD.Layout.Row = 6;
h.modeDD.Layout.Column = 2;
end


function [panel, h] = buildDICPanel(parent, p)
panel = uipanel(parent, 'Title', '2. DIC parameters', 'FontWeight', 'bold');
g = uigridlayout(panel, [4, 4]);
g.ColumnWidth = {200, 90, 200, 90};
g.RowHeight = repmat({28}, 1, 4);

uilabel(g, 'Text', 'Subset size (winsize, px):');
h.winsize = uispinner(g, 'Value', fallbackVal(p, 'winsize', 32), 'Limits', [4, 512], 'Step', 2);
h.winsize.Tooltip = 'Even integer. Physical DIC subset is (winsize+1) px wide. Typical 24-64 px.';

uilabel(g, 'Text', 'Step size (winstepsize, px):');
h.winstepsize = uispinner(g, 'Value', fallbackVal(p, 'winstepsize', 16), 'Limits', [1, 256], 'Step', 1);
h.winstepsize.Tooltip = 'Power of 2. Spacing between subset centers. Typical 8-32 px.';

uilabel(g, 'Text', 'Quadtree min size (winsizeMin, px):');
h.winsizeMin = uispinner(g, 'Value', fallbackVal(p, 'winsizeMin', 8), 'Limits', [1, 64], 'Step', 1);
h.winsizeMin.Tooltip = 'Smallest quadtree element near mask boundaries. Power of 2, <= winstepsize.';

uilabel(g, 'Text', 'Parallel workers (ClusterNo):');
h.clusterNo = uispinner(g, 'Value', fallbackVal(p, 'ClusterNo', 1), 'Limits', [1, 64], 'Step', 1);
h.clusterNo.Tooltip = '1 = serial. Higher values launch parpool. TemporalMatch auto-detects and enables 16 if left as 1.';

uilabel(g, 'Text', 'FFT search distance X (px):');
fft0 = safeFFTDist(p, 1);
h.fftDistX = uispinner(g, 'Value', fft0, 'Limits', [1, 500], 'Step', 5);
h.fftDistX.Tooltip = 'Initial-guess FFT search radius in X. Set larger than the maximum expected |disp u|.';

uilabel(g, 'Text', 'FFT search distance Y (px):');
fft1 = safeFFTDist(p, 2);
h.fftDistY = uispinner(g, 'Value', fft1, 'Limits', [1, 500], 'Step', 5);
h.fftDistY.Tooltip = 'Initial-guess FFT search radius in Y. Set larger than the maximum expected |disp v|.';

uilabel(g, 'Text', 'Use ALDIC global step:');
h.useGlobal = uicheckbox(g, 'Text', '', 'Value', fallbackVal(p, 'UseGlobal', true));
h.useGlobal.Tooltip = 'Checked: full ALDIC with ADMM (slower, more accurate). Unchecked: Local ICGN only (~2x faster, less smooth).';
end


function [panel, h] = buildStrainPanel(parent, p)
panel = uipanel(parent, 'Title', '3. Strain & smoothing', 'FontWeight', 'bold');
g = uigridlayout(panel, [2, 4]);
g.ColumnWidth = {200, 90, 200, 90};
g.RowHeight = repmat({28}, 1, 2);

uilabel(g, 'Text', 'Strain size (points):');
h.strainSize = uispinner(g, 'Value', fallbackVal(p, 'strain_size', 5), 'Limits', [2, 101], 'Step', 2);
h.strainSize.Tooltip = 'Side length of the local plane-fit window, in number of points. Larger = smoother strain, smaller = sharper features.';

uilabel(g, 'Text', 'Transform to specimen frame:');
h.transformDisp = uicheckbox(g, 'Text', '', 'Value', logical(fallbackVal(p, 'transformDisp', 0)));
h.transformDisp.Tooltip = 'If checked, prompts for 3 base points (O, X-axis, Y-axis) to rotate disp/strain into a specimen coordinate system.';

uilabel(g, 'Text', '2D smoothing count (ADMM):');
h.smooth2D = uispinner(g, 'Value', fallbackVal(p, 'Smooth2DTimes', 0), 'Limits', [0, 10], 'Step', 1);
h.smooth2D.Tooltip = 'Number of passes of inner-ADMM smoothing. 0 = disabled.';

uilabel(g, 'Text', '3D smoothing count (post-ALDIC):');
h.smooth3D = uispinner(g, 'Value', fallbackVal(p, 'Smooth3DTimes', 3), 'Limits', [0, 10], 'Step', 1);
h.smooth3D.Tooltip = 'Number of passes over FinalResult.Displacement (3D) before strain. 0 = raw disp, higher = smoother.';
end


function [panel, h] = buildOutputPanel(parent, p)
panel = uipanel(parent, 'Title', '4. Output', 'FontWeight', 'bold');
g = uigridlayout(panel, [3, 4]);
g.ColumnWidth = {200, 160, 200, 90};
g.RowHeight = repmat({28}, 1, 3);

uilabel(g, 'Text', 'Plot mode:');
h.plotMode = uidropdown(g, 'Items', {'show', 'save', 'none'}, ...
                       'Value', safeStr2(p, 'PlotMode', 'show'));
h.plotMode.Tooltip = 'show: live figures; save: render but do not display; none: no plots (future-proof; currently figures still render).';

uilabel(g, 'Text', 'Image to plot over:');
h.image2Plot = uidropdown(g, ...
    'Items', {'First frame only', 'Current (deformed) frame'}, ...
    'ItemsData', [0 1], ...
    'Value', fallbackVal(p, 'Image2PlotResults', 1));
h.image2Plot.Tooltip = 'Background image for the strain/displacement overlay.';

uilabel(g, 'Text', 'Figure format:');
h.figFormat = uidropdown(g, 'Items', {'png', 'jpg', 'pdf', 'eps', 'tif', 'fig'}, ...
                        'Value', safeStr2(p, 'MethodToSaveFig', 'png'));
h.figFormat.Tooltip = 'File type for saved figures.';

uilabel(g, 'Text', '');

uilabel(g, 'Text', 'Output folder:');
defaultOut = safeStr(p, 'outputFilePath');
if isempty(defaultOut)
    defaultOut = fullfile('.', 'results', char(datetime('now','Format','yyyyMMdd_HHmmss')));
end
h.outputPath = uieditfield(g, 'text', 'Value', defaultOut);
h.outputPath.Tooltip = 'Default: ./results/<timestamp>/. Created automatically if missing.';
browseOutBtn = uibutton(g, 'push', 'Text', 'Browse…', ...
    'ButtonPushedFcn', @(~,~) doPickOut());
browseOutBtn.Layout.Row = 3;
browseOutBtn.Layout.Column = 4;

    function doPickOut()
        folder = uigetdir(pwd, 'Pick output folder');
        if folder ~= 0
            h.outputPath.Value = folder;
        end
    end
end


% ================== Small utilities ==================

function addRow(g, rowIdx, labelText, browseFn)
uilabel(g, 'Text', labelText);
ed = uieditfield(g, 'text');
ed.Layout.Row = rowIdx; ed.Layout.Column = 2;
btn = uibutton(g, 'push', 'Text', 'Browse…', 'ButtonPushedFcn', ...
               @(~,~) setFromPicker(ed, browseFn));
btn.Layout.Row = rowIdx; btn.Layout.Column = 3;
end

function setFromPicker(edField, pickerFn)
val = pickerFn();
if ~isempty(val), edField.Value = val; end
end

function path = pickFolder(title)
folder = uigetdir(pwd, title);
if folder == 0, path = ''; else, path = folder; end
end

function path = pickFile(title)
[file, folder] = uigetfile({ ...
    '*.xml;*.mat;*.caldat;*.csv', 'Calibration files (*.xml,*.mat,*.caldat,*.csv)'; ...
    '*.*', 'All files (*.*)'}, title);
if isequal(file, 0), path = ''; else, path = fullfile(folder, file); end
end

function s = safeStr(p, field)
if isfield(p, field) && ~isempty(p.(field)) && (ischar(p.(field)) || isstring(p.(field)))
    s = char(p.(field));
else
    s = '';
end
end

function s = safeStr2(p, field, defaultVal)
if isfield(p, field) && ~isempty(p.(field))
    s = char(p.(field));
else
    s = defaultVal;
end
end

function v = fallbackVal(p, field, defaultVal)
if isfield(p, field) && ~isempty(p.(field))
    v = p.(field);
else
    v = defaultVal;
end
end

function d = safeFFTDist(p, idx)
if isfield(p, 'NewFFTSearchDistance') && ~isempty(p.NewFFTSearchDistance)
    arr = p.NewFFTSearchDistance;
    if length(arr) >= idx, d = arr(idx); else, d = arr(1); end
else
    d = 60;
end
end
