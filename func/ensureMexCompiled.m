function ensureMexCompiled(forceRebuild)
%ENSUREMEXCOMPILED  Compile ba_interp2_spline.cpp if missing or stale.
%
%   ensureMexCompiled()            % compile if missing or out of date
%   ensureMexCompiled(true)        % always recompile
%
%   Sets MW_MINGW64_LOC on Windows if a TDM-GCC install is in the default
%   location and the env var is not already defined. On any mex failure,
%   prints a platform-specific, actionable error message with concrete
%   install steps before raising.

if nargin < 1, forceRebuild = false; end

% On Windows, point MATLAB at a default-location TDM-GCC install.
if ispc && isempty(getenv('MW_MINGW64_LOC'))
    default_mingw = 'C:\TDM-GCC-64';
    if exist(default_mingw, 'dir')
        setenv('MW_MINGW64_LOC', default_mingw);
    end
end

mex_src = 'ba_interp2_spline.cpp';
mex_bin = ['ba_interp2_spline.', mexext];

if ~exist(mex_src, 'file')
    error('ensureMexCompiled:missingSource', ...
          'Cannot find %s. Run this from the 3D-Stereo-ALDIC project root.', mex_src);
end

needRebuild = forceRebuild || ~exist(mex_bin, 'file') || ...
              (dir(mex_src).datenum > dir(mex_bin).datenum);

if ~needRebuild
    fprintf('Mex binary up-to-date; skipping rebuild.\n');
    return;
end

fprintf('Compiling %s...\n', mex_src);
try
    mex('-O', mex_src);
    warning('off');
    fprintf('Mex compilation successful.\n');
catch ME
    printCompilerHelp(ME);
    error('ensureMexCompiled:failed', ...
          'MEX compile failed. Install a C++ compiler (see instructions above) and try again.');
end

end


function printCompilerHelp(ME)
%PRINTCOMPILERHELP  Emit a platform-specific guide for setting up mex.

fprintf(2, '\n');
fprintf(2, '============================================================\n');
fprintf(2, '  MEX compile failed\n');
fprintf(2, '============================================================\n');
fprintf(2, '  Original error:\n    %s\n\n', ME.message);

if ispc
    fprintf(2, '  Windows — set up a C++ compiler:\n');
    fprintf(2, '    Option 1 (MinGW, free, recommended):\n');
    fprintf(2, '      1) Download TDM-GCC-64:  https://jmeubank.github.io/tdm-gcc/\n');
    fprintf(2, '      2) Install to C:\\TDM-GCC-64 (default) or set the\n');
    fprintf(2, '         MW_MINGW64_LOC env var to your install path.\n');
    fprintf(2, '      3) In MATLAB: >> mex -setup C++ , pick MinGW.\n');
    fprintf(2, '\n');
    fprintf(2, '    Option 2 (Visual Studio Community, free):\n');
    fprintf(2, '      1) Install with the "Desktop development with C++" workload.\n');
    fprintf(2, '      2) In MATLAB: >> mex -setup C++ , pick Microsoft Visual C++.\n');
elseif ismac
    fprintf(2, '  macOS — set up a C++ compiler:\n');
    fprintf(2, '    1) In Terminal:  xcode-select --install\n');
    fprintf(2, '    2) In MATLAB:    >> mex -setup C++\n');
else
    fprintf(2, '  Linux — set up a C++ compiler:\n');
    fprintf(2, '    Ubuntu / Debian:  sudo apt install build-essential\n');
    fprintf(2, '    Fedora / RHEL:    sudo dnf install gcc-c++\n');
    fprintf(2, '    Then in MATLAB:   >> mex -setup C++\n');
end

fprintf(2, '\n');
fprintf(2, '  Full list of MATLAB-supported compilers:\n');
fprintf(2, '    https://www.mathworks.com/support/requirements/supported-compilers.html\n');
fprintf(2, '============================================================\n\n');

end
