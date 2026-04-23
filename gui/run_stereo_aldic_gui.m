function run_stereo_aldic_gui()
%RUN_STEREO_ALDIC_GUI  Tabbed GUI front-end for 3D-Stereo-ALDIC.
%
%   Opens a single uifigure with three tabs:
%     1. Configuration — all parameters, file pickers, Run button.
%     2. Log           — scrolling output from pipeline stages.
%     3. Visualize     — post-pipeline result explorer (enabled after Run):
%                        frame slider, variable dropdown, colorbar sliders,
%                        CSV export.
%
%   Algorithm / pipeline functions are unchanged; this is purely the
%   front-end.

%% ============ Environment setup ============
sessionStart = tic;

if isempty(getenv('MW_MINGW64_LOC'))
    default_mingw = 'C:\TDM-GCC-64';
    if exist(default_mingw, 'dir')
        setenv('MW_MINGW64_LOC', default_mingw);
    end
end

mex_src = 'ba_interp2_spline.cpp';
mex_bin = ['ba_interp2_spline.', mexext];
if ~exist(mex_bin, 'file') || (dir(mex_src).datenum > dir(mex_bin).datenum)
    try
        mex('-O', mex_src);
    catch ME
        warning('Mex compile failed: %s', ME.message);
    end
end

addpath('./examples', './func', './func_quadtree/rbfinterp/', ...
        './plotFiles/', './func_quadtree', './func_quadtree/refinement', ...
        './plotFiles/export_fig-d966721/', './gui', './gui/gui_helpers');

%% ============ State container ============
state.DICpara    = setDICparaDefaults();
state.results    = [];          % filled after Run
state.isRunning  = false;
state.logLines   = {};

%% ============ Main figure + tabs ============
fig = uifigure('Name', '3D-Stereo-ALDIC', 'Position', [120 80 1100 820]);
tabGroup = uitabgroup(fig, 'Position', [1 1 1100 820]);
tabConfig = uitab(tabGroup, 'Title', '  Configuration  ');
tabLog    = uitab(tabGroup, 'Title', '  Log  ');
tabViz    = uitab(tabGroup, 'Title', '  Visualize  ');

%% ============ Config tab widgets ============
hCfg = buildConfigTab(tabConfig, state.DICpara);

% Wire Config-tab callbacks
hCfg.runBtn.ButtonPushedFcn     = @(~,~) onRunClicked();
hCfg.viewCalibBtn.ButtonPushedFcn = @(~,~) onViewCalibration();
hCfg.drawROIBtn.ButtonPushedFcn = @(~,~) onDrawROI();

%% ============ Log tab widgets ============
hLog.area = uitextarea(tabLog, 'Position', [20 20 1060 760], ...
                      'Editable', 'off', 'FontName', 'Consolas', ...
                      'FontSize', 11, ...
                      'Value', {'Log ready. Hit Run on the Configuration tab to start.'});

%% ============ Visualize tab widgets ============
hViz = buildVisualizeTab(tabViz);
hViz.tabHandle = tabViz;
% Disable viz tab controls until results are ready
setVizEnabled(false);

% Wire viz callbacks
hViz.varDD.ValueChangedFcn    = @(~,~) updatePlot();
hViz.frameSlider.ValueChangedFcn = @(~,~) onFrameSlider();
hViz.cminSlider.ValueChangedFcn  = @(~,~) updatePlot();
hViz.cmaxSlider.ValueChangedFcn  = @(~,~) updatePlot();
hViz.autoScaleBtn.ButtonPushedFcn = @(~,~) autoScaleColors();
hViz.bgDD.ValueChangedFcn        = @(~,~) updatePlot();
hViz.exportBtn.ButtonPushedFcn   = @(~,~) onExportCSV();

% ============ Done. Show fig. ============

% ====================== Callbacks ======================

    function appendLog(msg)
        state.logLines{end+1} = msg;
        % Cap to last 500 lines to avoid runaway memory
        if length(state.logLines) > 500
            state.logLines = state.logLines(end-499:end);
        end
        hLog.area.Value = state.logLines;
        scroll(hLog.area, 'bottom');
        drawnow limitrate;
    end

    function onRunClicked()
        if state.isRunning
            return;
        end
        try
            state.DICpara = collectConfigValues();
            validateDICpara(state.DICpara);
            validateConfigPaths(state.DICpara);
        catch err
            hCfg.statusLabel.Text = ['ERROR: ' err.message];
            hCfg.statusLabel.FontColor = [0.8 0 0];
            return;
        end
        hCfg.statusLabel.Text = 'Running...';
        hCfg.statusLabel.FontColor = [0 0.4 0];
        hCfg.runBtn.Enable = 'off';
        state.isRunning = true;

        tabGroup.SelectedTab = tabLog;  % switch to log view
        state.logLines = {sprintf('=== Starting pipeline at %s ===', char(datetime('now')))};
        hLog.area.Value = state.logLines;
        drawnow;

        progressDlg = uiprogressdlg(fig, 'Title', 'Running pipeline', ...
                                    'Message', 'Starting...', 'Cancelable', 'on', 'Value', 0);

        function progressFn(frac, msg)
            if ~isvalid(progressDlg), error('User cancelled.'); end
            progressDlg.Value = frac;
            progressDlg.Message = msg;
        end

        try
            state.results = runPipelineCore(state.DICpara, @appendLog, @progressFn);
            appendLog(sprintf('=== Pipeline completed: %.1fs ===', toc(sessionStart)));
            hCfg.statusLabel.Text = sprintf('Done. Results in %s', state.DICpara.outputFilePath);
            hCfg.statusLabel.FontColor = [0 0.4 0];
            populateVisualize();
            setVizEnabled(true);
            tabGroup.SelectedTab = tabViz;
        catch err
            appendLog(sprintf('ERROR: %s', err.message));
            hCfg.statusLabel.Text = ['ERROR: ' err.message];
            hCfg.statusLabel.FontColor = [0.8 0 0];
        end

        if isvalid(progressDlg), close(progressDlg); end
        hCfg.runBtn.Enable = 'on';
        state.isRunning = false;
    end

    function DICpara = collectConfigValues()
        DICpara = state.DICpara;
        % Data
        DICpara.imgFolderLeft     = hCfg.imgLeftEdit.Value;
        DICpara.imgFolderRight    = hCfg.imgRightEdit.Value;
        DICpara.maskFolderLeft    = hCfg.maskLeftEdit.Value;
        DICpara.maskFolderRight   = hCfg.maskRightEdit.Value;
        DICpara.maskFolder        = DICpara.maskFolderLeft;  % legacy compat
        DICpara.calibrationFile   = hCfg.calibEdit.Value;
        DICpara.calibrationMethod = hCfg.calibFormatDD.Value;
        DICpara.DICIncOrNot       = hCfg.modeDD.Value;
        DICpara.manualROI         = state.manualROI;  % (may be empty)

        % DIC
        DICpara.winsize            = hCfg.winsize.Value;
        DICpara.winstepsize        = hCfg.winstepsize.Value;
        DICpara.winsizeMin         = hCfg.winsizeMin.Value;
        DICpara.ClusterNo          = hCfg.clusterNo.Value;
        DICpara.UseGlobal          = hCfg.useGlobal.Value;
        DICpara.StereoSearchDistance   = hCfg.stereoSearch.Value;   % scalar -> square
        DICpara.TemporalSearchDistance = hCfg.temporalSearch.Value;

        % Strain
        DICpara.strain_size    = hCfg.strainSize.Value;
        DICpara.Smooth2DTimes  = hCfg.smooth2D.Value;
        DICpara.Smooth3DTimes  = hCfg.smooth3D.Value;
        DICpara.transformDisp  = double(hCfg.transformDisp.Value);

        % Output
        DICpara.PlotMode          = hCfg.plotMode.Value;
        DICpara.Image2PlotResults = hCfg.image2Plot.Value;
        DICpara.MethodToSaveFig   = hCfg.figFormat.Value;
        DICpara.outputFilePath    = hCfg.outputPath.Value;
    end

    function onViewCalibration()
        fpath = hCfg.calibEdit.Value;
        if isempty(fpath) || ~exist(fpath, 'file')
            uialert(fig, 'Pick a calibration file first.', 'Calibration', 'Icon', 'warning');
            return;
        end
        method = hCfg.calibFormatDD.Value;
        [folder, name, ext] = fileparts(fpath);
        try
            info = struct();
            switch method
                case 0, info = cameraParamsFormatConvertFromMatlabCV;
                case 1, info.cameraParams = cameraParamsFormatConvertFromMatchID(folder, [name ext]);
                case 2, info.cameraParams = cameraParamsFormatConvertFromMMC(folder, [name ext]);
                case 3, info.cameraParams = cameraParamsFormatConvertFromDICe(folder, [name ext]);
                case 4, info.cameraParams = cameraParamsFormatConvertFromOpenCorrFormat(folder, [name ext]);
            end
            showCalibrationDialog(fig, info.cameraParams);
        catch err
            uialert(fig, sprintf('Could not parse calibration:\n%s', err.message), ...
                    'Calibration error', 'Icon', 'error');
        end
    end

    function onDrawROI()
        if hCfg.modeDD.Value == 1
            uialert(fig, ['Manual ROI drawing is only supported in Accumulative mode. ' ...
                          'For Incremental mode, please use the mask-folder inputs.'], ...
                    'ROI drawing', 'Icon', 'info');
            return;
        end
        leftFolder = hCfg.imgLeftEdit.Value;
        if isempty(leftFolder) || ~exist(leftFolder, 'dir')
            uialert(fig, 'Pick a valid Left images folder first.', 'ROI', 'Icon', 'warning');
            return;
        end
        [~, firstImg] = pickFirstImage(leftFolder);
        if isempty(firstImg), return; end
        state.manualROI = drawROIOnImage(firstImg);
        if ~isempty(state.manualROI)
            hCfg.statusLabel.Text = 'Manual ROI captured (overrides mask bounding box).';
            hCfg.statusLabel.FontColor = [0 0.4 0];
        end
    end

    % ---- Visualize tab callbacks ----
    function populateVisualize()
        r = state.results;
        if isempty(r), return; end
        nFrames = size(r.FinalResult.Coordinates, 1);

        hViz.frameSlider.Limits = [2 nFrames];
        hViz.frameSlider.Value  = 2;
        hViz.frameSlider.MajorTicks = 2:max(1,floor((nFrames-2)/5)):nFrames;
        hViz.frameLabel.Text = sprintf('Frame: 2 / %d', nFrames);

        % Reset colorbar sliders (auto-scale on first plot)
        hViz.cminSlider.Value = 0;
        hViz.cmaxSlider.Value = 1;

        updatePlot();
    end

    function updatePlot()
        if isempty(state.results), return; end
        r = state.results;
        frameIdx = round(hViz.frameSlider.Value);
        varName  = hViz.varDD.Value;
        bgChoice = hViz.bgDD.Value;

        hViz.frameLabel.Text = sprintf('Frame: %d / %d', frameIdx, size(r.FinalResult.Coordinates, 1));

        % Pick data vector per variable
        [vals, displayName] = getVariableValues(r, frameIdx, varName);
        if isempty(vals)
            cla(hViz.ax); title(hViz.ax, 'No data'); return;
        end

        coords = r.RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM;
        if size(coords, 1) ~= length(vals)
            % Incremental mode: strain is on FirstFEM, but coords in the
            % FirstFEM ~ always aligned. Still fall back to vals size.
        end

        % Background image
        bgImg = [];
        try
            if bgChoice == 1  % current frame
                fname = fullfile(r.fileNameLeft{2,frameIdx}, r.fileNameLeft{1,frameIdx});
            else
                fname = fullfile(r.fileNameLeft{2,1}, r.fileNameLeft{1,1});
            end
            raw = imread(fname);
            if size(raw,3) >= 3, raw = rgb2gray(raw(:,:,1:3)); end
            bgImg = raw;
        catch
            bgImg = [];
        end

        % Color scale
        finite = vals(isfinite(vals));
        if isempty(finite)
            vmin = 0; vmax = 1;
        else
            [vmin, vmax] = deal(min(finite), max(finite));
        end
        % Apply user slider positions (0..1 of vmin..vmax range)
        range = vmax - vmin;
        cmin = vmin + hViz.cminSlider.Value * range;
        cmax = vmin + hViz.cmaxSlider.Value * range;
        if cmax <= cmin, cmax = cmin + eps(max(abs(cmin), 1)); end

        cla(hViz.ax);
        if ~isempty(bgImg)
            imshow(bgImg, 'Parent', hViz.ax);
            hold(hViz.ax, 'on');
        end
        scatter(hViz.ax, coords(:,1), coords(:,2), 6, vals, 'filled', ...
                'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.85);
        hold(hViz.ax, 'off');
        colormap(hViz.ax, 'turbo');
        cb = colorbar(hViz.ax);
        cb.Label.String = displayName;
        clim(hViz.ax, [cmin cmax]);
        axis(hViz.ax, 'image');
        set(hViz.ax, 'YDir', 'reverse');
        title(hViz.ax, sprintf('%s — frame %d', displayName, frameIdx), ...
              'Interpreter', 'none');
    end

    function onFrameSlider()
        hViz.frameSlider.Value = round(hViz.frameSlider.Value);
        updatePlot();
    end

    function autoScaleColors()
        hViz.cminSlider.Value = 0;
        hViz.cmaxSlider.Value = 1;
        updatePlot();
    end

    function onExportCSV()
        if isempty(state.results)
            uialert(fig, 'Run the pipeline first.', 'Export', 'Icon', 'warning');
            return;
        end
        [file, folder] = uiputfile({'*.csv','CSV (*.csv)'}, 'Export current view', ...
                                   fullfile(state.DICpara.outputFilePath, ...
                                            sprintf('%s_frame%d.csv', hViz.varDD.Value, round(hViz.frameSlider.Value))));
        if isequal(file, 0), return; end
        r = state.results;
        frameIdx = round(hViz.frameSlider.Value);
        varName = hViz.varDD.Value;
        [vals, ~] = getVariableValues(r, frameIdx, varName);
        coords = r.RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM;
        n = min(length(vals), size(coords,1));
        T = table(coords(1:n,1), coords(1:n,2), vals(1:n), ...
                  'VariableNames', {'X_image', 'Y_image', varName});
        writetable(T, fullfile(folder, file));
        appendLog(sprintf('Exported %d rows to %s', n, fullfile(folder, file)));
    end

    function setVizEnabled(flag)
        children = [hViz.varDD, hViz.frameSlider, hViz.cminSlider, ...
                    hViz.cmaxSlider, hViz.bgDD, hViz.autoScaleBtn, hViz.exportBtn];
        for k = 1:length(children)
            if flag, children(k).Enable = 'on'; else, children(k).Enable = 'off'; end
        end
    end

end  % main function


% ============================================================
% Helper: build Config tab
% ============================================================
function h = buildConfigTab(parent, p)
mainGrid = uigridlayout(parent, [6 1]);
mainGrid.RowHeight = {180, 160, 120, 130, 'fit', 40};
mainGrid.RowSpacing = 8; mainGrid.Padding = [14 14 14 14];

% --- Data ---
panel1 = uipanel(mainGrid, 'Title', '1. Data', 'FontWeight', 'bold');
g1 = uigridlayout(panel1, [6 3]);
g1.ColumnWidth = {220, '1x', 90};
g1.RowHeight = repmat({28}, 1, 6);

[h.imgLeftEdit,  ~] = addPathRow(g1, 1, 'Left images folder:',      p, 'imgFolderLeft',      @pickFolder);
[h.imgRightEdit, ~] = addPathRow(g1, 2, 'Right images folder:',     p, 'imgFolderRight',     @pickFolder);
[h.maskLeftEdit, ~] = addPathRow(g1, 3, 'Left mask folder:',        p, 'maskFolderLeft',     @pickFolder);
[h.maskRightEdit,~] = addPathRow(g1, 4, 'Right mask folder:',       p, 'maskFolderRight',    @pickFolder);
[h.calibEdit,    ~] = addPathRow(g1, 5, 'Calibration file:',        p, 'calibrationFile',    @pickCalibFile);

uilabel(g1, 'Text', 'Calibration format + Mode:');
formatGrid = uigridlayout(g1, [1 2]);
formatGrid.ColumnSpacing = 6; formatGrid.Padding = [0 0 0 0];
formatGrid.Layout.Row = 6; formatGrid.Layout.Column = 2;
h.calibFormatDD = uidropdown(formatGrid, ...
    'Items', {'MATLAB CV (0)','MatchID (1)','MCC (2)','DICe XML (3)','OpenCorr (4)'}, ...
    'ItemsData', [0 1 2 3 4], ...
    'Value', valOr(p, 'calibrationMethod', 3));
h.modeDD = uidropdown(formatGrid, ...
    'Items', {'Accumulative (small disp)','Incremental (large disp)'}, ...
    'ItemsData', [0 1], ...
    'Value', valOr(p, 'DICIncOrNot', 0));
calExtraBtns = uigridlayout(g1, [1 2]);
calExtraBtns.Layout.Row = 6; calExtraBtns.Layout.Column = 3;
calExtraBtns.Padding = [0 0 0 0];
calExtraBtns.ColumnSpacing = 3;
h.viewCalibBtn = uibutton(calExtraBtns, 'push', 'Text', 'View calib');
h.drawROIBtn   = uibutton(calExtraBtns, 'push', 'Text', 'Draw ROI');

% --- DIC ---
panel2 = uipanel(mainGrid, 'Title', '2. DIC parameters', 'FontWeight', 'bold');
g2 = uigridlayout(panel2, [4 4]);
g2.ColumnWidth = {220, 90, 220, 90};
g2.RowHeight = repmat({28}, 1, 4);

uilabel(g2, 'Text', 'Subset size (winsize, px):');
h.winsize = uispinner(g2, 'Value', valOr(p,'winsize',32), 'Limits', [4 512], 'Step', 2);
h.winsize.Tooltip = 'Even integer. Physical subset = (winsize+1) px.';

uilabel(g2, 'Text', 'Step size (winstepsize, px):');
h.winstepsize = uispinner(g2, 'Value', valOr(p,'winstepsize',16), 'Limits', [1 256], 'Step', 1);
h.winstepsize.Tooltip = 'Power of 2. Spacing between subsets.';

uilabel(g2, 'Text', 'Quadtree min size (winsizeMin):');
h.winsizeMin = uispinner(g2, 'Value', valOr(p,'winsizeMin',8), 'Limits', [1 64], 'Step', 1);

uilabel(g2, 'Text', 'Parallel workers:');
h.clusterNo = uispinner(g2, 'Value', valOr(p,'ClusterNo',1), 'Limits', [1 64], 'Step', 1);

uilabel(g2, 'Text', 'Stereo FFT search radius (px):');
h.stereoSearch = uispinner(g2, 'Value', valOrFirst(p, 'StereoSearchDistance', 80), ...
                           'Limits', [1 500], 'Step', 5);
h.stereoSearch.Tooltip = 'Square search radius for the stereo match (frame 1 only). Set larger than max expected disparity.';

uilabel(g2, 'Text', 'Temporal FFT search radius (px):');
h.temporalSearch = uispinner(g2, 'Value', valOrFirst(p, 'TemporalSearchDistance', 20), ...
                             'Limits', [1 500], 'Step', 5);
h.temporalSearch.Tooltip = 'Square search radius for frame-to-frame match. Set larger than max inter-frame motion.';

uilabel(g2, 'Text', 'Use ALDIC global step:');
h.useGlobal = uicheckbox(g2, 'Text', '', 'Value', valOr(p,'UseGlobal',true));

% --- Strain ---
panel3 = uipanel(mainGrid, 'Title', '3. Strain & smoothing', 'FontWeight', 'bold');
g3 = uigridlayout(panel3, [2 4]);
g3.ColumnWidth = {220, 90, 220, 90}; g3.RowHeight = {28 28};
uilabel(g3, 'Text', 'Strain size (points):');
h.strainSize = uispinner(g3, 'Value', valOr(p,'strain_size',5), 'Limits', [2 101], 'Step', 2);

uilabel(g3, 'Text', 'Transform to specimen frame:');
h.transformDisp = uicheckbox(g3, 'Text', '', 'Value', logical(valOr(p,'transformDisp',0)));

uilabel(g3, 'Text', '2D smoothing count:');
h.smooth2D = uispinner(g3, 'Value', valOr(p,'Smooth2DTimes',0), 'Limits', [0 10]);

uilabel(g3, 'Text', '3D smoothing count:');
h.smooth3D = uispinner(g3, 'Value', valOr(p,'Smooth3DTimes',3), 'Limits', [0 10]);

% --- Output ---
panel4 = uipanel(mainGrid, 'Title', '4. Output', 'FontWeight', 'bold');
g4 = uigridlayout(panel4, [3 4]);
g4.ColumnWidth = {220, 160, 220, 90}; g4.RowHeight = {28 28 28};
uilabel(g4, 'Text', 'Plot mode (during pipeline):');
h.plotMode = uidropdown(g4, 'Items', {'none','save','show'}, ...
                        'Value', strOr(p,'PlotMode','none'));
h.plotMode.Tooltip = 'none = skip per-frame plots (fastest, recommended; use Visualize tab to inspect results).';

uilabel(g4, 'Text', 'Image to plot over:');
h.image2Plot = uidropdown(g4, 'Items', {'First frame only','Current (deformed) frame'}, ...
                          'ItemsData', [0 1], 'Value', valOr(p,'Image2PlotResults',1));

uilabel(g4, 'Text', 'Figure format:');
h.figFormat = uidropdown(g4, 'Items', {'png','jpg','pdf','eps','tif','fig'}, ...
                         'Value', strOr(p,'MethodToSaveFig','png'));
uilabel(g4, 'Text', '');

uilabel(g4, 'Text', 'Output folder:');
defaultOut = strOr(p,'outputFilePath', fullfile('.', 'results', char(datetime('now','Format','yyyyMMdd_HHmmss'))));
h.outputPath = uieditfield(g4, 'text', 'Value', defaultOut);
outBrowseBtn = uibutton(g4, 'push', 'Text', 'Browse…', ...
                       'ButtonPushedFcn', @(~,~) (function_handle_wrap(@pickFolder, h.outputPath)));
outBrowseBtn.Layout.Row = 3; outBrowseBtn.Layout.Column = 4;

% --- Status ---
h.statusLabel = uilabel(mainGrid, 'Text', '', 'WordWrap', 'on');

% --- Run button ---
runBar = uigridlayout(mainGrid, [1 3]);
runBar.ColumnWidth = {'1x', 140, 140};
uilabel(runBar, 'Text', 'Ready. Click Run ▶ to start.');
cancelBtn = uibutton(runBar, 'push', 'Text', 'Cancel GUI');
cancelBtn.ButtonPushedFcn = @(~,~) closereq_fig();
h.runBtn = uibutton(runBar, 'push', 'Text', 'Run ▶', ...
                   'BackgroundColor', [0.2 0.6 0.2], 'FontColor', 'white', 'FontWeight', 'bold');
end


% ============================================================
% Helper: build Visualize tab
% ============================================================
function h = buildVisualizeTab(parent)
g = uigridlayout(parent, [2 2]);
g.RowHeight    = {'1x', 120};
g.ColumnWidth  = {'1x', 280};
g.Padding      = [10 10 10 10];

% Main axes
h.ax = uiaxes(g);
h.ax.Layout.Row = 1; h.ax.Layout.Column = [1 2];

% Controls on bottom row
controlPanel = uipanel(g, 'Title', 'View controls');
controlPanel.Layout.Row = 2; controlPanel.Layout.Column = [1 2];
cg = uigridlayout(controlPanel, [3 6]);
cg.ColumnWidth = {90, 140, 90, '1x', 90, 120};
cg.RowHeight   = {28, 28, 28};

% Variable + background
uilabel(cg, 'Text', 'Variable:');
h.varDD = uidropdown(cg, 'Items', ...
    {'U','V','W','|U|','epsilon_xx','epsilon_yy','epsilon_xy','epsilon_1','epsilon_2','maxshear','vonMises'}, ...
    'Value', 'epsilon_xx');
uilabel(cg, 'Text', 'Background:');
h.bgDD = uidropdown(cg, 'Items', {'Current frame','First frame'}, ...
                   'ItemsData', [1 0], 'Value', 1);
h.bgDD.Layout.Row = 1; h.bgDD.Layout.Column = 4;
h.autoScaleBtn = uibutton(cg, 'push', 'Text', 'Auto color');
h.autoScaleBtn.Layout.Row = 1; h.autoScaleBtn.Layout.Column = 5;
h.exportBtn = uibutton(cg, 'push', 'Text', 'Export CSV…');
h.exportBtn.Layout.Row = 1; h.exportBtn.Layout.Column = 6;

% Frame slider
uilabel(cg, 'Text', 'Frame:');
h.frameSlider = uislider(cg, 'Limits', [2 3], 'Value', 2);
h.frameSlider.Layout.Row = 2; h.frameSlider.Layout.Column = [2 5];
h.frameLabel = uilabel(cg, 'Text', 'Frame: 2 / 2');
h.frameLabel.Layout.Row = 2; h.frameLabel.Layout.Column = 6;

% Colorbar sliders
uilabel(cg, 'Text', 'Color range:');
h.cminSlider = uislider(cg, 'Limits', [0 1], 'Value', 0);
h.cminSlider.Layout.Row = 3; h.cminSlider.Layout.Column = 2;
uilabel(cg, 'Text', 'min ↔ max');
h.cmaxSlider = uislider(cg, 'Limits', [0 1], 'Value', 1);
h.cmaxSlider.Layout.Row = 3; h.cmaxSlider.Layout.Column = 4;
end


% ============================================================
% Misc helpers
% ============================================================
function [vals, displayName] = getVariableValues(r, frameIdx, varName)
displayName = varName;
vals = [];
try
    switch varName
        case 'U', vals = r.FinalResult.Displacement{frameIdx, 1};
        case 'V', vals = r.FinalResult.Displacement{frameIdx, 2};
        case 'W', vals = r.FinalResult.Displacement{frameIdx, 3};
        case '|U|'
            u = r.FinalResult.Displacement{frameIdx,1};
            v = r.FinalResult.Displacement{frameIdx,2};
            w = r.FinalResult.Displacement{frameIdx,3};
            vals = sqrt(u.^2 + v.^2 + w.^2);
        case {'epsilon_xx','epsilon_yy','epsilon_xy','epsilon_1','epsilon_2','maxshear','vonMises'}
            s = r.strainPerFrame{frameIdx};
            if isempty(s), return; end
            switch varName
                case 'epsilon_xx', vals = s.strain_exx;
                case 'epsilon_yy', vals = s.strain_eyy;
                case 'epsilon_xy', vals = s.strain_exy;
                case 'epsilon_1',  vals = s.strain_principal_max;
                case 'epsilon_2',  vals = s.strain_principal_min;
                case 'maxshear',   vals = s.strain_maxshear;
                case 'vonMises',   vals = s.strain_vonMises;
            end
    end
catch
    vals = [];
end
end


function [ed, btn] = addPathRow(g, rowIdx, labelText, p, field, pickFn)
uilabel(g, 'Text', labelText);
ed = uieditfield(g, 'text', 'Value', strOr(p, field, ''));
btn = uibutton(g, 'push', 'Text', 'Browse…', 'ButtonPushedFcn', @(~,~) doPick());
    function doPick()
        v = pickFn();
        if ~isempty(v), ed.Value = v; end
    end
end

function path = pickFolder()
folder = uigetdir(pwd);
if isequal(folder, 0), path = ''; else, path = folder; end
end

function path = pickCalibFile()
[f, folder] = uigetfile({'*.xml;*.mat;*.caldat;*.csv','Calibration'; '*.*','All'});
if isequal(f, 0), path = ''; else, path = fullfile(folder, f); end
end

function function_handle_wrap(pickFn, field)
v = pickFn();
if ~isempty(v), field.Value = v; end
end

function v = valOr(p, field, def)
if isfield(p, field) && ~isempty(p.(field)), v = p.(field); else, v = def; end
end

function v = valOrFirst(p, field, def)
if isfield(p, field) && ~isempty(p.(field))
    a = p.(field); v = a(1);
else
    v = def;
end
end

function v = strOr(p, field, def)
if isfield(p, field) && ~isempty(p.(field)), v = char(p.(field)); else, v = def; end
end

function validateConfigPaths(p)
if isempty(p.imgFolderLeft)  || ~exist(p.imgFolderLeft,'dir'),  error('Left images folder missing.'); end
if isempty(p.imgFolderRight) || ~exist(p.imgFolderRight,'dir'), error('Right images folder missing.'); end
if isempty(p.maskFolderLeft) || ~exist(p.maskFolderLeft,'dir'), error('Left mask folder missing.'); end
if isempty(p.maskFolderRight)|| ~exist(p.maskFolderRight,'dir'),error('Right mask folder missing.'); end
if isempty(p.calibrationFile)|| ~exist(p.calibrationFile,'file'),error('Calibration file missing.'); end
end

function showCalibrationDialog(~, cp)
d = uifigure('Name','Calibration details','Position',[200 200 520 500]);
ta = uitextarea(d,'Position',[10 10 500 480],'Editable','off','FontName','Consolas','FontSize',11);
lines = {};
lines = append_mat(lines, 'Left camera K',  tryGet(cp, 'cameraParamsLeft.K'));
lines = append_mat(lines, 'Right camera K', tryGet(cp, 'cameraParamsRight.K'));
lines = append_mat(lines, 'Rotation matrix (L -> R)', tryGet(cp, 'rotationMatrix'));
lines = append_mat(lines, 'Translation vector (mm)',  tryGet(cp, 'translationVector'));
lines = append_mat(lines, 'Left distortion',  tryGet(cp, 'cameraParamsLeft.RadialDistortion'));
lines = append_mat(lines, 'Right distortion', tryGet(cp, 'cameraParamsRight.RadialDistortion'));
ta.Value = lines;
end

function lines = append_mat(lines, label, val)
if isempty(val), return; end
lines{end+1} = '';
lines{end+1} = ['=== ' label ' ==='];
if isnumeric(val) && ~isscalar(val)
    lines{end+1} = mat2str(val, 6);
else
    lines{end+1} = sprintf('%.5g ', val);
end
end

function v = tryGet(obj, dotPath)
v = [];
parts = strsplit(dotPath, '.');
cur = obj;
for k = 1:length(parts)
    if isfield(cur, parts{k}) || isprop(cur, parts{k})
        cur = cur.(parts{k});
    else
        return;
    end
end
v = cur;
end

function [fullPath, img] = pickFirstImage(folder)
exts = {'jpg','jpeg','tif','tiff','bmp','png','jp2'};
entries = [];
for k=1:length(exts)
    entries = [entries; dir(fullfile(folder, ['*.' exts{k}]))];  %#ok<AGROW>
end
if isempty(entries), fullPath=''; img=[]; return; end
[~, idx] = sort({entries.name}); entries = entries(idx);
fullPath = fullfile(entries(1).folder, entries(1).name);
img = imread(fullPath);
if size(img,3)>=3, img = rgb2gray(img(:,:,1:3)); end
end

function roi = drawROIOnImage(img)
roi = [];
f = figure('Name','Draw ROI polygon, double-click to finish');
imshow(img); hold on;
try
    poly = drawpolygon('LineWidth',2,'Color','y');
    wait(poly);
    roi = poly.Position;
    close(f);
catch
    try, close(f); catch, end
end
end

function closereq_fig()
delete(gcbf);
end
