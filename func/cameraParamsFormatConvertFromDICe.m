function cameraParams_matlab = cameraParamsFormatConvertFromDICe(CalibrationFilepath, CalibrationFile)
% This function is used to convert calibration format from DICe (cal.xml file)

cameraParams_Nodes = xmlread(fullfile(CalibrationFilepath, CalibrationFile));

paramLists = cameraParams_Nodes.getElementsByTagName('ParameterList');

% Initialize partial variable to store the extracted value
LeftFs = 0;
RightFs = 0;

% Traverse all ParameterList nodes to find "CAMERA 0"
for i = 0:paramLists.getLength-1
    paramList = paramLists.item(i);

    % Left camera
    if strcmp(char(paramList.getAttribute('name')), 'CAMERA 0')
        % Find all Parameter nodes in "CAMERA 0"
        parameters = paramList.getElementsByTagName('Parameter');

        % Traverse the Parameter nodes to find "CX"...
        for k = 0:parameters.getLength-1
            parameterNode = parameters.item(k);
            paramName = char(parameterNode.getAttribute('name'));
            paramValue = str2double(char(parameterNode.getAttribute('value')));

            % Extract value based on parameter name
            switch paramName
                case 'CX'
                    LeftCx = paramValue;
                case 'CY'
                    LeftCy = paramValue;
                case 'FX'
                    LeftFx = paramValue;
                case 'FY'
                    LeftFy = paramValue;
                case 'FS'
                    LeftFs = paramValue;
                case 'K1'
                    LeftK1 = paramValue;
                case 'K2'
                    LeftK2 = paramValue;
                case 'K3'
                    LeftK3 = paramValue;
                case 'P1'
                    LeftP1 = paramValue;
                case 'P2'
                    LeftP2 = paramValue;
            end
        end
    end

    % Right camera
    if strcmp(char(paramList.getAttribute('name')), 'CAMERA 1')
        % Find all Parameter nodes in "CAMERA 1"
        parameters = paramList.getElementsByTagName('Parameter');

        % Traverse the Parameter nodes to find "CX"...
        for k = 0:parameters.getLength-1
            parameterNode = parameters.item(k);
            paramName = char(parameterNode.getAttribute('name'));

            % Extract value based on parameter name
            switch paramName
                case 'CX'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightCx = paramValue;
                case 'CY'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightCy = paramValue;
                case 'FX'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightFx = paramValue;
                case 'FY'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightFy = paramValue;
                case 'FS'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightFs = paramValue;
                case 'K1'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightK1 = paramValue;
                case 'K2'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightK2 = paramValue;
                case 'K3'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightK3 = paramValue;
                case 'P1'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightP1 = paramValue;
                case 'P2'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    RightP2 = paramValue;
                case 'TX'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    Tx = paramValue;
                case 'TY'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    Ty = paramValue;
                case 'TZ'
                    paramValue = str2double(char(parameterNode.getAttribute('value')));
                    Tz = paramValue;
                case 'ROW 0'
                    paramValue = char(parameterNode.getAttribute('value'));
                    cleanStr = strrep(paramValue, '{', '[');
                    cleanStr = strrep(cleanStr, '}', ']');
                    cleanStr = strrep(cleanStr, ' ', '');
                    ROW_0 = str2num(cleanStr);
                case 'ROW 1'
                    paramValue = char(parameterNode.getAttribute('value'));
                    cleanStr = strrep(paramValue, '{', '[');
                    cleanStr = strrep(cleanStr, '}', ']');
                    cleanStr = strrep(cleanStr, ' ', '');
                    ROW_1 = str2num(cleanStr);
                case 'ROW 2'
                    paramValue = char(parameterNode.getAttribute('value'));
                    cleanStr = strrep(paramValue, '{', '[');
                    cleanStr = strrep(cleanStr, '}', ']');
                    cleanStr = strrep(cleanStr, ' ', '');
                    ROW_2 = str2num(cleanStr);
            end
        end
    end
end

% Conversion
cameraParams_matlab.cameraParamsLeft.K = [LeftFx LeftFs LeftCx; 0 LeftFy LeftCy ; 0 0 1];
cameraParams_matlab.cameraParamsRight.K = [RightFx RightFs RightCx; 0 RightFy RightCy ; 0 0 1];
cameraParams_matlab.rotationMatrix = [ROW_0;ROW_1;ROW_2];
cameraParams_matlab.translationVector = [Tx,Ty,Tz];

try
cameraParams_matlab.cameraParamsLeft.RadialDistortion = [LeftK1,LeftK2,LeftK3];
cameraParams_matlab.cameraParamsRight.RadialDistortion = [RightK1,RightK2,RightK3];
catch
cameraParams_matlab.cameraParamsLeft.RadialDistortion = [0,0,0];
cameraParams_matlab.cameraParamsRight.RadialDistortion = [0,0,0];
end

try 
cameraParams_matlab.cameraParamsLeft.TangentialDistortion = [LeftP1,LeftP2];
cameraParams_matlab.cameraParamsRight.TangentialDistortion = [RightP1,RightP2];
catch
cameraParams_matlab.cameraParamsLeft.TangentialDistortion = [0,0];
cameraParams_matlab.cameraParamsRight.TangentialDistortion = [0,0];
end

