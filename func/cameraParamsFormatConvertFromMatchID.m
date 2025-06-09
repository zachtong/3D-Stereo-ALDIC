function cameraParams_matlab = cameraParamsFormatConvertFromMatchID(CalibrationFilepath, CalibrationFile)
% This function is used to convert calibration format in matchid
% to that in matlab.

cameraParams_matchid = importdata(fullfile(CalibrationFilepath, CalibrationFile));

% Assign
LeftFx = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Fx")));
LeftFy = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Fy")));
LeftFs = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Fs")));
LeftCx = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Cx")));
LeftCy = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Cy")));


RightFx = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Fx")));
RightFy = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Fy")));
RightFs = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Fs")));
RightCx = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Cx")));
RightCy = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Cy")));

try  % Adapt different calib. output format
    LeftK1 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa1")),2);
    LeftK2 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa2")),2);
    LeftK3 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa3")),2);
    LeftP1 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam0_P1")),2);
    LeftP2 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam0_P2")),2);
    RightK1 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa1")),2);
    RightK2 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa2")),2);
    RightK3 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa3")),2);
    RightP1 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam1_P1")),2);
    RightP2 = cameraParams_matchid.textdata(find(contains(cameraParams_matchid.textdata,"Cam1_P2")),2);
catch
    LeftK1 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa 1")));
    LeftK2 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa 2")));
    LeftK3 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_Kappa 3")));
    LeftP1 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_P1")));
    LeftP2 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam0_P2")));
    RightK1 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa 1")));
    RightK2 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa 2")));
    RightK3 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_Kappa 3")));
    RightP1 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_P1")));
    RightP2 = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Cam1_P2")));
end

Tx = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Tx")));
Ty = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Ty")));
Tz = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Tz")));

Theta = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Theta")));
Phi = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Phi")));
Psi = cameraParams_matchid.data(find(contains(cameraParams_matchid.textdata,"Psi")));

% Conversion
cameraParams_matlab.cameraParamsLeft.K = [LeftFx LeftFs LeftCx; 0 LeftFy LeftCy ; 0 0 1];
cameraParams_matlab.cameraParamsRight.K = [RightFx RightFs RightCx; 0 RightFy RightCy ; 0 0 1];
cameraParams_matlab.rotationMatrix = rotation_matrix(Theta, Phi, Psi);
cameraParams_matlab.translationVector = [Tx,Ty,Tz];

try % Adapt different calib. output format
    cameraParams_matlab.cameraParamsLeft.RadialDistortion = [str2double(LeftK1{1}),str2double(LeftK2{1}),str2double(LeftK3{1})];
    cameraParams_matlab.cameraParamsLeft.TangentialDistortion = [str2double(LeftP1),str2double(LeftP2)];
    cameraParams_matlab.cameraParamsRight.RadialDistortion = [str2double(RightK1{1}),str2double(RightK2{1}),str2double(RightK3{1})];
    cameraParams_matlab.cameraParamsRight.TangentialDistortion = [str2double(RightP1),str2double(RightP2)];
catch
    cameraParams_matlab.cameraParamsLeft.RadialDistortion = [LeftK1,LeftK2,LeftK3];
    cameraParams_matlab.cameraParamsLeft.TangentialDistortion = [LeftP1,LeftP2];
    cameraParams_matlab.cameraParamsRight.RadialDistortion = [RightK1,RightK2,RightK3];
    cameraParams_matlab.cameraParamsRight.TangentialDistortion = [RightP1,RightP2];
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