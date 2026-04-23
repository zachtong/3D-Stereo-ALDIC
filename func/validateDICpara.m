function validateDICpara(DICpara)
%VALIDATEDICPARA  Validate DICpara fields, error on misconfiguration.
%
%   validateDICpara(DICpara)
%
%   Checks types, ranges, and mutual consistency of DICpara fields.
%   Raises an MException with a clear message if any field is invalid.
%
%   Mirrors pyALDIC's al_dic.core.config.validate_dicpara.
%
%   Typical usage:
%       DICpara = setDICparaDefaults();
%       DICpara.winsize     = 32;
%       DICpara.winstepsize = 16;
%       validateDICpara(DICpara);
%
%   See also setDICparaDefaults, main_3D_ALDIC.

isPow2 = @(v) isnumeric(v) && isscalar(v) && v > 0 && bitand(int32(v), int32(v) - 1) == 0;

% ---- required-for-run fields (skipped if empty; caller sets later) ----
if ~isempty(DICpara.winsize)
    if ~isnumeric(DICpara.winsize) || ~isscalar(DICpara.winsize) || ...
            DICpara.winsize <= 0 || mod(DICpara.winsize, 2) ~= 0
        error('validateDICpara:badWinsize', ...
            'winsize must be a positive even integer, got %g', DICpara.winsize);
    end
end

if ~isempty(DICpara.winstepsize)
    if ~isPow2(DICpara.winstepsize)
        error('validateDICpara:badWinstep', ...
            'winstepsize must be a positive power of 2, got %g', DICpara.winstepsize);
    end
    if ~isPow2(DICpara.winsizeMin)
        error('validateDICpara:badWinsizeMin', ...
            'winsizeMin must be a positive power of 2, got %g', DICpara.winsizeMin);
    end
    if DICpara.winsizeMin > DICpara.winstepsize
        error('validateDICpara:badWinsizeMin', ...
            'winsizeMin (%d) must be <= winstepsize (%d)', ...
            DICpara.winsizeMin, DICpara.winstepsize);
    end
end

% ---- enumerations ----
if ~ismember(DICpara.DICIncOrNot, [0 1])
    error('validateDICpara:badIncOrNot', ...
        'DICIncOrNot must be 0 (acc) or 1 (inc), got %g', DICpara.DICIncOrNot);
end

if ~ismember(DICpara.StrainType, 0:3)
    error('validateDICpara:badStrainType', ...
        'StrainType must be 0-3 (currently only 0 = Green-Lagrange supported)');
end

if ~ischar(DICpara.PlotMode) && ~isstring(DICpara.PlotMode)
    error('validateDICpara:badPlotMode', 'PlotMode must be a string');
end
if ~ismember(char(DICpara.PlotMode), {'none', 'save', 'show'})
    error('validateDICpara:badPlotMode', ...
        'PlotMode must be ''none'', ''save'', or ''show'', got ''%s''', ...
        char(DICpara.PlotMode));
end

if ~ismember(DICpara.InitFFTSearchMethod, [0 1 2 3])
    error('validateDICpara:badFFTMethod', ...
        'InitFFTSearchMethod must be 0-3, got %g', DICpara.InitFFTSearchMethod);
end

% ---- smoothing counts ----
if ~isscalar(DICpara.Smooth2DTimes) || DICpara.Smooth2DTimes < 0 || ...
        mod(DICpara.Smooth2DTimes, 1) ~= 0
    error('validateDICpara:badSmooth2D', ...
        'Smooth2DTimes must be a non-negative integer, got %g', DICpara.Smooth2DTimes);
end
if ~isscalar(DICpara.Smooth3DTimes) || DICpara.Smooth3DTimes < 0 || ...
        mod(DICpara.Smooth3DTimes, 1) ~= 0
    error('validateDICpara:badSmooth3D', ...
        'Smooth3DTimes must be a non-negative integer, got %g', DICpara.Smooth3DTimes);
end

% ---- cluster / misc ----
if ~isscalar(DICpara.ClusterNo) || DICpara.ClusterNo < 1
    error('validateDICpara:badClusterNo', ...
        'ClusterNo must be >= 1, got %g', DICpara.ClusterNo);
end

if ~islogical(DICpara.UseGlobal) && ~ismember(DICpara.UseGlobal, [0 1])
    error('validateDICpara:badUseGlobal', ...
        'UseGlobal must be logical or 0/1');
end

if ~isempty(DICpara.NewFFTSearchDistance)
    d = DICpara.NewFFTSearchDistance;
    if ~isnumeric(d) || ~any(numel(d) == [1 2 4]) || any(d <= 0)
        error('validateDICpara:badFFTDist', ...
            'NewFFTSearchDistance must be [], scalar, [x y], or [L R U D], all positive');
    end
end

if ~isempty(DICpara.progressCallback) && ~isa(DICpara.progressCallback, 'function_handle')
    error('validateDICpara:badCallback', ...
        'progressCallback must be a function handle @(frac, msg) or empty');
end

end
