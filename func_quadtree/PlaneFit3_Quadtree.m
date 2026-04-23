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

try
    R = DICpara.RforStrainCal;  % If DICpara.RforStrainCal exists, that means we already define a new coor sys
catch
    Strain_Coordinate = 'Local';
end

Coordinate3D_reference = [Coordinate3D_reference{1},Coordinate3D_reference{2},Coordinate3D_reference{3}];

 
%% Strain calculation
dilatedI = logical(DICpara.ImgRefMask);
cc = bwconncomp(dilatedI,8);
[row1,~] = find(round(AllPoints2D(:,1))>1);
[row2,~] = find(round(AllPoints2D(:,1))<DICpara.ImgSize(1));
[row3,~] = find(round(AllPoints2D(:,2))>1);
[row4,~] = find(round(AllPoints2D(:,2))<DICpara.ImgSize(2));
row1234 = intersect(intersect(intersect(row1,row2),row3),row4);


% indPxAll: linear indices of all node pixels in the reference image
indPxAll = sub2ind( DICpara.ImgSize, round(AllPoints2D(row1234,1)), round(AllPoints2D(row1234,2)) );

% One 'stats' entry per connected mask region; strain is computed
% independently inside each region.
stats = regionprops(cc,'Area','PixelList');
for connAreaNum = 1:length(stats)
    % indPxtempi: linear indices of all pixels in this connected region
    indPxtempi = sub2ind( DICpara.ImgSize, stats(connAreaNum).PixelList(:,2), stats(connAreaNum).PixelList(:,1) );

    % Pick the nodes that fall inside this region
    Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);

    % Per-node displacement + pixel coordinates for this region
    tempDispUVW = Displacement(LiaList,:);
    tempCoor = AllPoints2D(LiaList,:);
    tempCoordinate3D_reference = Coordinate3D_reference(LiaList,:);
    % Strain is computed at every node (strain winstep = 1). Knn_number is
    % an upper bound on the neighborhood we search; scale with VSG so we
    % always find enough points inside the window.
    Knn_number = 2 * ceil(strain_length/DICpara.winsizeMin)^2;
    if Knn_number < 9, Knn_number = 9; end  % minimum 9 nodes for stable plane fit
    [neighborInd,distance] = knnsearch(tempCoor,tempCoor,'K',Knn_number,'Distance','chebychev');

    % Keep only neighbors inside the VSG (Chebyshev radius = 0.5 * strain_length).
    % Boundary nodes of a region may have fewer neighbors; STAQ accepts this
    % tradeoff for higher edge resolution.
    % Note: max(.., [],2) on a logical row returns the index of the first
    % 'true' (first OUTSIDE neighbor). minus-1 yields the count of neighbors
    % INSIDE the VSG. If all neighbors are inside, max returns 1 and we end
    % up with 0 — a known edge case (would only trigger for tiny datasets).
    OutsideOrNot = distance > 0.5 * strain_length;
    [~,numInsideVSG] = max(OutsideOrNot, [], 2);
    numInsideVSG = numInsideVSG - 1;

    if all(numInsideVSG < 9), disp('Please increase VSG size!'); end

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

