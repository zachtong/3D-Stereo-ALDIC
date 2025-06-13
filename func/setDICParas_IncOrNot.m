function DICpara = setDICParas_IncOrNot(numImages)
if numImages > 2
    % ==============================================
    % DIC initial guess
    % NewFFTSearch = funParaInput('NewFFTSearch'); % Use last frame as init guess or not

    % ==============================================
    % IncrementalOrNot  #
    IncrementalOrNot = funParaInput('IncrementalOrNot');

    try
        switch IncrementalOrNot
            case 0 % acc. default
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
                NewFFTSearch = 0;
            case 1 % inc.
                ImgSeqIncUnit = 1;
                ImgSeqIncROIUpdateOrNot = 0; 
                NewFFTSearch = 1;
            otherwise
                ImgSeqIncUnit = numImages+1;
                ImgSeqIncROIUpdateOrNot = 1;
        end
    catch
        ImgSeqIncUnit = numImages+1;
        ImgSeqIncROIUpdateOrNot = 1;
    end

    % ================================
else % Only two frames
    ImgSeqIncUnit = numImages+1;
    ImgSeqIncROIUpdateOrNot = 1;
    IncrementalOrNot = 0; % Revised by Zach
end

% Double check
NewFFTSearch =1; %tbd
DICpara.NewFFTSearch = NewFFTSearch;

DICpara.DICIncOrNot = IncrementalOrNot; % Revised by Zach
DICpara.ImgSeqIncUnit = ImgSeqIncUnit;
DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot;
end