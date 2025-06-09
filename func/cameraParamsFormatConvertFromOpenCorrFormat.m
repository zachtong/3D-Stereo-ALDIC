function cameraParams_OpenCorr = cameraParamsFormatConvertFromOpenCorrFormat(CalibrationFilepath, CalibrationFile)
% This function is used to convert calibration format from OpenCorr (.csv file)
data = readtable(fullfile(CalibrationFilepath, CalibrationFile));

cameraParams_OpenCorr = struct();

for i = 1:height(data)
    variableName = data{i, 1}{1};
    value_left = data{i, 2};
    value_right = data{i, 3};

    % Extract value based on parameter name
    switch variableName
        case 'Cx'
            LeftCx = value_left; RightCx = value_right;
        case 'Cy'
            LeftCy = value_left; RightCy = value_right;
        case 'Fx'
            LeftFx = value_left; RightFx = value_right;
        case 'Fy'
            LeftFy = value_left; RightFy = value_right;
        case 'Fs'
            LeftFs = value_left; RightFs = value_right;
        case 'K1'
            LeftK1 = value_left; RightK1 = value_right;
        case 'K2'
            LeftK2 = value_left; RightK2 = value_right;
        case 'K3'
            LeftK3 = value_left; RightK3 = value_right;
        case 'K4'
            LeftK4 = value_left; RightK4 = value_right;
        case 'K5'
            LeftK5 = value_left; RightK5 = value_right;
        case 'K6'
            LeftK6 = value_left; RightK6 = value_right;
        case 'P1'
            LeftP1 = value_left; RightP1 = value_right;
        case 'P2'
            LeftP2 = value_left; RightP2 = value_right;
        case 'Tx'
            RightTx = value_right;
        case 'Ty'
            RightTy = value_right;
        case 'Tz'
            RightTz = value_right;
        case 'Rx'
            RightRx = value_right;
        case 'Ry'
            RightRy = value_right;
        case 'Rz'
            RightRz = value_right;
    end
end

% Conversion
cameraParams_OpenCorr.cameraParamsLeft.K = [LeftFx LeftFs LeftCx; 0 LeftFy LeftCy ; 0 0 1];
cameraParams_OpenCorr.cameraParamsRight.K = [RightFx RightFs RightCx; 0 RightFy RightCy ; 0 0 1];
cameraParams_OpenCorr.rotationMatrix = rotation_matrix(RightRx, RightRy, RightRz);
cameraParams_OpenCorr.translationVector = [RightTx,RightTy,RightTz];

try
    cameraParams_OpenCorr.cameraParamsLeft.RadialDistortion = [LeftK1,LeftK2,LeftK3];
    cameraParams_OpenCorr.cameraParamsRight.RadialDistortion = [RightK1,RightK2,RightK3];
catch
    cameraParams_OpenCorr.cameraParamsLeft.RadialDistortion = [0,0,0];
    cameraParams_OpenCorr.cameraParamsRight.RadialDistortion = [0,0,0];
end

try
    cameraParams_OpenCorr.cameraParamsLeft.TangentialDistortion = [LeftP1,LeftP2];
    cameraParams_OpenCorr.cameraParamsRight.TangentialDistortion = [RightP1,RightP2];
catch
    cameraParams_OpenCorr.cameraParamsLeft.TangentialDistortion = [0,0];
    cameraParams_OpenCorr.cameraParamsRight.TangentialDistortion = [0,0];
end
end


function R = rotation_matrix(Theta, Phi, Psi) % Support Unit[degree]
Rz = [cosd(Psi), -sind(Psi), 0;
    sind(Psi), cosd(Psi), 0;
    0, 0, 1];

Ry = [cosd(Phi), 0, sind(Phi);
    0, 1, 0;
    -sind(Phi), 0, cosd(Phi)];

Rx = [1, 0, 0;
    0, cosd(Theta), -sind(Theta);
    0, sind(Theta), cosd(Theta)];

R = Rz * Ry * Rx;
end
