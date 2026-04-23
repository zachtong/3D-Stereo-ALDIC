function DICpara = setDICParas_IncOrNot(numImages)
%SETDICPARAS_INCORNOT  Ask user for accumulative vs. incremental mode.
%
%   Returns a partial DICpara struct containing DICIncOrNot,
%   ImgSeqIncUnit, ImgSeqIncROIUpdateOrNot, and NewFFTSearch. Other
%   fields are filled by setDICParas_STAQ.

if numImages > 2
    IncrementalOrNot = funParaInput('IncrementalOrNot');
    switch IncrementalOrNot
        case 0  % accumulative (default)
            ImgSeqIncUnit = numImages + 1;
            ImgSeqIncROIUpdateOrNot = 1;
        case 1  % incremental
            ImgSeqIncUnit = 1;
            ImgSeqIncROIUpdateOrNot = 0;
        otherwise
            ImgSeqIncUnit = numImages + 1;
            ImgSeqIncROIUpdateOrNot = 1;
    end
else  % only two frames
    ImgSeqIncUnit = numImages + 1;
    ImgSeqIncROIUpdateOrNot = 1;
    IncrementalOrNot = 0;
end

DICpara.DICIncOrNot             = IncrementalOrNot;
DICpara.ImgSeqIncUnit           = ImgSeqIncUnit;
DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot;
DICpara.NewFFTSearch            = 1;  % always re-run FFT initial guess
end
