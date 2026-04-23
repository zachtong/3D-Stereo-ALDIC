function DICpara = setDICparaDefaults()
%SETDICPARADEFAULTS  Return a DICpara struct with all fields initialized.
%
%   DICpara = setDICparaDefaults()
%
%   Central source of truth for every DICpara field used anywhere in the
%   3D-Stereo-ALDIC codebase.  Fields that MUST be set by the caller are
%   initialized to [] (empty).  Fields with sensible defaults are filled.
%
%   Pattern (mirrors pyALDIC's al_dic.core.config.dicpara_default):
%       DICpara = setDICparaDefaults();
%       DICpara.winsize     = 32;        % user override
%       DICpara.winstepsize = 16;
%       validateDICpara(DICpara);        % catch misconfig
%
%   See also validateDICpara, main_3D_ALDIC.

DICpara = struct();

% ---------------------------------------------------------------------
% Core DIC algorithm parameters (user MUST set winsize / winstepsize)
% ---------------------------------------------------------------------
DICpara.winsize              = [];      % subset half-width (code convention, physical size = winsize+1)
DICpara.winstepsize          = [];      % subset sampling step, power of 2
DICpara.winsizeMin           = 8;       % min element size in quadtree mesh
DICpara.DICIncOrNot          = 0;       % 0 = accumulative, 1 = incremental
DICpara.NewFFTSearch         = 1;       % re-run FFT initial guess each frame
DICpara.UseGlobal            = true;    % run ADMM global step (Subpb2)
DICpara.Subpb2FDOrFEM        = 1;       % 1 = FEM, 2 = finite difference (currently only 1 supported in quadtree path)
DICpara.GaussPtOrder         = 2;
DICpara.ClusterNo            = 1;       % parpool size; 1 = serial

% ---------------------------------------------------------------------
% FFT initial-guess search region (set to [] to force interactive prompt)
% ---------------------------------------------------------------------
DICpara.NewFFTSearchDistance   = [];    % e.g. [60, 60] px for stereo disparity
DICpara.fixSearchDistanceOrNot = 0;     % 0 = reuse preset for all frames, 1 = re-prompt per frame

% ---------------------------------------------------------------------
% Image / ROI
% ---------------------------------------------------------------------
DICpara.ImgRefMask           = [];      % set by caller from mask image
DICpara.ImgSize              = [];      % set by caller from image size
DICpara.gridxyROIRange       = struct('gridx', [], 'gridy', []);
DICpara.InitFFTSearchMethod  = 1;       % 1 = whole-field, 2 = seeded, 3 = zero, 0 = multigrid
DICpara.LoadImgMethod        = 1;

% ---------------------------------------------------------------------
% Strain calculation
% ---------------------------------------------------------------------
DICpara.strain_size          = [];      % side length of local plane-fit window (pts); prompts if empty
DICpara.transformDisp        = 0;       % 0 = keep camera coords, 1 = transform to specimen via GetRTMatrix
DICpara.um2px                = 1;       % unit conversion factor
DICpara.StrainType           = 0;       % only 0 (Green-Lagrange) is implemented

% ---------------------------------------------------------------------
% Smoothing (independently controllable counts)
% ---------------------------------------------------------------------
DICpara.Smooth2DTimes        = 0;       % inside ADMM (2D disp/strain); 0 = disabled
DICpara.Smooth3DTimes        = 3;       % post-ALDIC (3D disp on FinalResult); 3 = current default
DICpara.DispSmoothness       = 0;       % curvature regularizer weight (in funSmoothDispQuadtree)
DICpara.StrainSmoothness     = 0;
DICpara.DispFilterSize       = 0;
DICpara.DispFilterStd        = 0;
DICpara.StrainFilterSize     = 0;
DICpara.StrainFilterStd      = 0;

% ---------------------------------------------------------------------
% Plotting / output
% ---------------------------------------------------------------------
DICpara.PlotMode             = 'show';  % 'none' | 'save' | 'show'
DICpara.showImgOrNot         = 0;       % debug figures inside algorithm functions
DICpara.plots_disp_to_generate   = {'u', 'v', 'w', 'magnitude'};
DICpara.plots_strain_to_generate = {'exx', 'eyy', 'exy', 'e1', 'e2'};
DICpara.Image2PlotResults    = 1;       % 0 = first image only, 1 = second+
DICpara.MethodToSaveFig      = 'png';
DICpara.OrigDICImgTransparency = 0.8;
DICpara.outputFilePath       = '';      % filled with timestamped path at runtime if empty

% ---------------------------------------------------------------------
% MEX build control
% ---------------------------------------------------------------------
DICpara.forceMexRebuild      = false;

% ---------------------------------------------------------------------
% Optional progress callback: @(frac, msg) fprintf(...)
% ---------------------------------------------------------------------
DICpara.progressCallback     = [];

end
