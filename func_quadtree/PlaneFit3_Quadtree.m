function [DICpara,coefficients,voidIndex] = PlaneFit3_Quadtree( strain_size, DICpara, AllPoints2D ,...
    Coordinate3D_reference, Coordinate3D_current, Displacement,FullImageName_first, ...
    imageNum ,Strain_Coordinate)
%% This function is used for local plane fitting for the coefficients of U,V,W (eg: U = ax + by + cz + d)
% Note: Variable Coordinate_ref is used for Specified Strain Coordinate.
% User can assign 3 points as O,X,Y points to generate O-XYZ coordinate.

manuallySelectOrNot = 0;

strain_length = (strain_size-1) * DICpara.winstepsize + 1;
coefficients = cell(size(Coordinate3D_reference{1},1),1); % Ux Uy Uz


% For Lagrange config.
Displacement = [Displacement{1},Displacement{2},Displacement{3}];

% For Euler config.
% Coordinate3D_current = [Coordinate3D_current{1},Coordinate3D_current{2},Coordinate3D_current{3}];
% Displacement = [Displacement{1},Displacement{2},Displacement{3}];

%% Preparation for Specific Coordinate: Convert to a particular Coor. sys.
% if imageNum == 26
%     switch Strain_Coordinate
%         % Only for approximate plane objects
%         case 'Specific Coordinate'
%             if manuallySelectOrNot == 0
%                 imshow(FullImageName_first);
%                 title('Select three points as O,X,Y base points','fontweight','normal','fontsize',16);
% 
%                 Point_O = zeros(1,2); Point_X = zeros(1,2); Point_Y = zeros(1,2);
%                 [Point_O(1), Point_O(2)] = ginput(1);
%                 fprintf( 'Point O: (%4.3f,%4.3f) \n' ,Point_O(1), Point_O(2))
%                 [Point_X(1), Point_X(2)] = ginput(1);
%                 fprintf( 'Point X: (%4.3f,%4.3f) \n' ,Point_X(1), Point_X(2))
%                 [Point_Y(1), Point_Y(2)] = ginput(1);
%                 fprintf( 'Point Y: (%4.3f,%4.3f) \n' ,Point_Y(1), Point_Y(2))
% 
%             else
%                 % % Users manual input three points
%                 % % Ch2 task1
%                 % Point_O = [866,2368];
%                 % Point_X = [1255,2355];
%                 % Point_Y = [844,784];
% 
%                 % Point_O = [980,1700];
%                 % Point_X = [1110,1700];
%                 % Point_Y = [980,1400];
% 
%                 % Set new X- Y- Z- axis by choosing 3 points manually
%                 % Left-upper corner is (0,0). first element represents the number of row.
%                 Point_O = [496,541];
%                 Point_X = [900,541];
%                 Point_Y = [496,158];
% 
% 
%             end
% 
%             Base_Points_2D = [Point_O;Point_X;Point_Y];
%             R = GetRTMatrix(AllPoints2D,Base_Points_2D,Coordinate3D_reference);
%             DICpara.RforStrainCal = R;
%     end
% else % For imageNum > 2, no need to repeat redundant calculations
%     try
%         R = DICpara.RforStrainCal;
%     catch
%         Strain_Coordinate = 'Local';
%     end
% end

try
    R = DICpara.RforStrainCal;  % If DICpara.RforStrainCal exists, that means we already define a new coor sys
catch
    Strain_Coordinate = 'Local';
end

Coordinate3D_reference = [Coordinate3D_reference{1},Coordinate3D_reference{2},Coordinate3D_reference{3}];

 
%% Strain calculation
dilatedI = logical(DICpara.ImgRefMask); % figure, imshow(dilatedI);
cc = bwconncomp(dilatedI,8);
[row1,~] = find(round(AllPoints2D(:,1))>1);
[row2,~] = find(round(AllPoints2D(:,1))<DICpara.ImgSize(1));
[row3,~] = find(round(AllPoints2D(:,2))>1);
[row4,~] = find(round(AllPoints2D(:,2))<DICpara.ImgSize(2));
row1234 = intersect(intersect(intersect(row1,row2),row3),row4);
% figure, plot(round(AllPoints2D(row1234,1)),round(AllPoints2D(row1234,2)),'.');


% indPxAll 是所有nodes的pixel对应的index
indPxAll = sub2ind( DICpara.ImgSize, round(AllPoints2D(row1234,1)), round(AllPoints2D(row1234,2)) );

% 每个Stat代表一个连通域，每个连通域分开计算内部strain
stats = regionprops(cc,'Area','PixelList');
for connAreaNum = 1:length(stats)
    % figure, plot(stats(connAreaNum).PixelList(:,1), stats(connAreaNum).PixelList(:,2),'.');

    % indPxtempi 是该连通域内所有的pixel对应的index
    indPxtempi = sub2ind( DICpara.ImgSize, stats(connAreaNum).PixelList(:,2), stats(connAreaNum).PixelList(:,1) );

    % 先找 indPxAll 中有哪些nodes是在该连通域内
    Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);

    % 位于该连通域内的 nodes 位移和像素坐标
    tempDispUVW = Displacement(LiaList,:);
    tempCoor = AllPoints2D(LiaList,:);  
    tempCoordinate3D_reference = Coordinate3D_reference(LiaList,:);
    % 确定每个 nodes 的 Knn_number 个最近邻 nodes
    % 目前每个node都计算应变，相当于strain winstep = 1
    Knn_number = 2 * ceil(strain_length/DICpara.winsizeMin)^2;
    if Knn_number < 9; Knn_number = 9; else; end % At least search 9 nodes
    [neighborInd,distance] = knnsearch(tempCoor,tempCoor,'K',Knn_number,'Distance','chebychev');

    % 在这 Knn_number 个最近邻 nodes 中，筛选出位于VSG内的 nodes
    OutsideOrNot = distance > 0.5 * strain_length;
    [~,numInsideVSG] = max(OutsideOrNot, [], 2);
    numInsideVSG = numInsideVSG - 1; % Actually only # of (numInsideVSG-1) inside

    % 是否去掉位于连通域边缘的 nodes？ （判断依据： numInsideVSG < strain_size^2）
    % 暂时不去，STAQ的目的就是为了增加边缘区域的空间分辨率，有足够的点就可以算
    % 唯一需要注意的地方是，对于这些边缘的点，他们的 VSG size 不等于 Overall VSG size

    % 如果所有的VSG内都少于 9 个 numInsideVSG
    if all(numInsideVSG < 9); disp('Please increase VSG size!'); end

    % Calculate strain for each VSG
    for i = 1:size(numInsideVSG,1) 
        tempIndex = neighborInd(i,1:numInsideVSG(i));
        % Get a,b,c,d by Least squares method (U/V/W = ax + by + cz +d)
        % [subDisp] = [subMatrix] .* [Coef]
        subMatrix = [tempCoordinate3D_reference(tempIndex,:)]; % Lagrange strain: reference_coor
        subDisp = tempDispUVW(tempIndex,:); % ~Dimension: (2*R+1)^2,1
        switch Strain_Coordinate
            case 'Camera0'
                % LSM
                Camera0_Corr_Coef = [subMatrix,ones(numInsideVSG(i), 1)]\subDisp;
                % Save data
                coefficients_temp = Camera0_Corr_Coef(1:3,:);
            case 'Local' %% Step 2 - Option 1: Convert dudx to Local Coor.
                % What was previously fitted was the four-dimensional
                % surface of UVW, but now what needs to be re-projected
                % is the three-dimensional surface of the actual shape. Z = ax + by + c

                subXY = [tempCoordinate3D_reference(tempIndex,1:2),ones(numInsideVSG(i), 1)];
                subZ = tempCoordinate3D_reference(tempIndex,3);
                Local_Corr_Coef = subXY\subZ;

                % plane_fitting_error_mean(i,:) = mean(abs(subXY * Local_Corr_Coef - subZ));
                % plane_fitting_error_max(i,:) = max(abs(subXY * Local_Corr_Coef - subZ));

                % Rotation matrix
                direction_z = [Local_Corr_Coef(1); Local_Corr_Coef(2) ;-1]; direction_z = direction_z/norm(direction_z);
                direction_x = [1;0;0] - [1;0;0]'*direction_z*direction_z; direction_x = direction_x/norm(direction_x);
                direction_y = cross(direction_z,direction_x);
                R = [direction_x,direction_y,direction_z];

                % R11{i,1} = R(1,1);
                % R22{i,1} = R(2,2);
                % R33{i,1} = R(3,3);

                % VSG center
                % Coor_VSG_center = subMatrix((square_size^2+1)/2,:);
                % Coor_VSG_center = Coordinate(tempIndex((square_size^2+1)/2),:); % 新的局部平面建立用的是deformed_coor
                % Local_subMatrix = ( (subMatrix - Coor_VSG_center) * R);

                Local_subMatrix = ( (subMatrix ) * R);

                % 策略1：直接用 U = aX + bY + cZ + d 来拟合，然后再将所有的c都置为0
                % Local_subDisp =  subDisp * R;
                % Local_Corr_Coef = [Local_subMatrix,ones(numInsideVSG(i), 1)]\Local_subDisp;

                % 策略2：直接用 U = aX + bY + d 来拟合，即视作所有的Z都是0，这也是一种策略
                % 但实际上，Z相比于X和Y本身就非常的小，几乎不影响拟合的参数
                Local_subDisp =  subDisp * R;
                Local_Corr_Coef = [Local_subMatrix(:,1:2),ones(numInsideVSG(i), 1)]\Local_subDisp;

                % Save data
                %coefficients( idx, :) = Local_Corr_Coef(1:3)';
                coefficients_temp = [Local_Corr_Coef(1:2,:); 0 0 0];
            case 'Specific Coordinate' %% Step 2 - Option 2: Convert to a particular Coor.

                Specific_subMatrix = ((subMatrix)*R );
                Specific_subDisp =  subDisp * R;
                Specific_Corr_Coef = [Specific_subMatrix(:,1:2),ones(numInsideVSG(i), 1)]\Specific_subDisp;
                % Save data
                % coefficients( idx, :) = Local_Corr_Coef(1:3)';
                coefficients_temp = [Specific_Corr_Coef(1:2,:); 0 0 0];
        end
            
        coefficients{LiaList(i),1} = coefficients_temp;
    end
end

% Fill void
voidIndex = cellfun('isempty', coefficients);
coefficients(voidIndex,1) = {zeros(3)};

end

