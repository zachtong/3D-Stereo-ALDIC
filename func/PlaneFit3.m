function [coefficients] = PlaneFit3(M, N, Rad, DICpara, Pixel_Coor_ref, Coordinate_ref, Coordinate, Displacement, Strain_Coordinate,task)
%% This function is used for local plane fitting for the coefficients of U,V,W (eg: U = ax + by + cz + d)
% Note: Variable Coordinate_ref is used for Specified Strain Coordinate.
% User can assign 3 points as O,X,Y points to generate O-XYZ coordinate.

square_size = 2*Rad +1;
indexMap = reshape(1:M*N, [M, N])';
Coordinate_ref = [Coordinate_ref{1},Coordinate_ref{2},Coordinate_ref{3}];
% Coordinate = [Coordinate{1},Coordinate{2},Coordinate{3}];

%% Preparation for Option 2: Convert to a particular Coor
switch Strain_Coordinate
    case 'Specific Coordinate'
        % Set new X- Y- Z- axis by choosing 3 points manually
        % Left-upper corner is (0,0). first element represents the number of row.

        %% ---------------------comment temporally
        % Ch2 task1
        % Point_O = [866,2368];
        % Point_X = [1255,2355];
        % Point_Y = [844,784];
        %
        % % Point_O = [980,1700];
        % % Point_X = [1110,1700];
        % % Point_Y = [980,1400];
        %
        %
        %
        % % Method 1: rbf interp (Global interpolation, high memory requirement)
        % % interpolation OXY corr in World_corr
        % rbfInterpX = rbfcreate((Pixel_Coor_ref)',(Coordinate_ref(:,1))');
        % rbfInterpY = rbfcreate((Pixel_Coor_ref)',(Coordinate_ref(:,2))');
        % rbfInterpZ = rbfcreate((Pixel_Coor_ref)',(Coordinate_ref(:,3))');
        %
        % % interpolating
        % Base_Points(1,:) = rbfinterp([Point_O',Point_X',Point_Y'],rbfInterpX);
        % Base_Points(2,:) = rbfinterp([Point_O',Point_X',Point_Y'],rbfInterpY);
        % Base_Points(3,:) = rbfinterp([Point_O',Point_X',Point_Y'],rbfInterpZ);
        % Base_Points = Base_Points';
        %%----------------------------------------------------------


        % % Method 2: scatteredInterpolant (faster, suitable for fewer calculation points)
        % Base_Points_pixel_coor = [Point_O;Point_X;Point_Y];
        % interp_x = scatteredInterpolant(Pixel_Coor_ref(:,1),Pixel_Coor_ref(:,2),Coordinate_ref(:,1),'natural');
        % interp_y = scatteredInterpolant(Pixel_Coor_ref(:,1),Pixel_Coor_ref(:,2),Coordinate_ref(:,2),'natural');
        % interp_z = scatteredInterpolant(Pixel_Coor_ref(:,1),Pixel_Coor_ref(:,2),Coordinate_ref(:,3),'natural');
        %
        % Base_Points(1,:) = interp_x(Base_Points_pixel_coor(:,1),Base_Points_pixel_coor(:,2));
        % Base_Points(2,:) = interp_y(Base_Points_pixel_coor(:,1),Base_Points_pixel_coor(:,2));
        % Base_Points(3,:) = interp_z(Base_Points_pixel_coor(:,1),Base_Points_pixel_coor(:,2));
        % Base_Points = Base_Points';
        %


        %% TASK 1&2 3 ---------------------------------------------
        if task ==1 || task ==2 || task ==3
            Base_Points = [-1.816023675315637,53.722142503543765,5.484058571239904e+02;
                17.915165549301953,53.358670285121110,5.512737314326490e+02;
                -2.945763928151216,-26.101729933259065,5.455083542765537e+02];
        elseif task ==4
            % TASK 4 ---------------------------------------------
            Base_Points = [-3.46486743240349	20.0354070015234	546.944886127772;
                14.3489433673661	19.8152026735436	549.472757910099;
                -3.54325271234572	-54.0755319894924	544.470353529227];
            % %%---------------------------------
        end

        %%
        new_x_axis = [Base_Points(2,:)-Base_Points(1,:)]'; new_x_axis = new_x_axis / norm(new_x_axis);
        new_y_axis = [Base_Points(3,:)-Base_Points(1,:)]'; new_y_axis = new_y_axis / norm(new_y_axis);
        new_z_axis = cross(new_x_axis,new_y_axis);

        % Rotation matrix between 2 coor. sys.
        R = [new_x_axis,new_y_axis,new_z_axis];
        disp(['TASK: ']);
        disp(task);


end


%% For each strain calculation window
if DICpara.ClusterNo <= 1 %% No-parfor
    idx = 1;
    for i = 1:N
        for j = 1:M

            if ~((i<1+Rad) || (i>N-Rad) || (j<1+Rad) || (j>M-Rad))
                %% Step 1 Get the dudx in Camera 1's Coor.
                % U/V/W = ax + by + cz +d ; Get a,b,c,d by Least squares method
                % Get the index of points in VSG
                Index = getSquareIDs(M, N, i, j, Rad, indexMap);

                %  [subDisp] = [subMatrix] .* [Coef]
                subMatrix = [Coordinate_ref(Index,:)]; % Lagrange strain: reference_coor

                subDisp = [Displacement{1,1},Displacement{1,2},Displacement{1,3}]; % Dimension: (2*R+1)^2,1
                subDisp = subDisp(Index,:);
                switch Strain_Coordinate
                    case 'Camera0'
                        % LSM
                        Camera0_Corr_Coef = [subMatrix,ones(square_size^2, 1)]\subDisp;
                        % Save data
                        coefficients{idx, :} = Camera0_Corr_Coef(1:3);
                        idx = idx+1;

                    case 'Local' %% Step 2 - Option 1: Convert dudx to Local Coor.
                        % What was previously fitted was the four-dimensional
                        % surface of UVW, but now what needs to be re-projected
                        % is the three-dimensional surface of the actual shape. Z = ax + b y + c

                        subXY = [Coordinate_ref(Index,1:2),ones(square_size^2, 1)];
                        subZ = Coordinate_ref(Index,3);
                        Local_Corr_Coef = subXY\subZ;

                        % Rotation matrix
                        direction_z = [Local_Corr_Coef(1); Local_Corr_Coef(2) ;-1]; direction_z = direction_z/norm(direction_z);
                        direction_x = [1;0;0] - [1;0;0]'*direction_z*direction_z; direction_x = direction_x/norm(direction_x);
                        direction_y = cross(direction_z,direction_x);
                        R = [direction_x,direction_y,direction_z];

                        % VSG center
                        % Coor_VSG_center = subMatrix((square_size^2+1)/2,:);
                        % Coor_VSG_center = Coordinate(Index((square_size^2+1)/2),:); % 新的局部平面建立用的是deformed_coor
                        % Local_subMatrix = ( (subMatrix - Coor_VSG_center) * R);

                        Local_subMatrix = ( (subMatrix ) * R);

                        % 策略1：直接用 U = aX + bY + cZ + d 来拟合，然后再将所有的c都置为0
                        % Local_subDisp =  subDisp * R;
                        % Local_Corr_Coef = [Local_subMatrix,ones(square_size^2, 1)]\Local_subDisp;


                        % 策略2：直接用 U = aX + bY + d 来拟合，即视作所有的Z都是0，这也是一种策略
                        % 但实际上，Z相比于X和Y本身就非常的小，几乎不影响拟合的参数
                        Local_subDisp =  subDisp * R;
                        Local_Corr_Coef = [Local_subMatrix(:,1:2),ones(square_size^2, 1)]\Local_subDisp;


                        % Save data
                        %coefficients( idx, :) = Local_Corr_Coef(1:3)';
                        coefficients{idx,1} = [Local_Corr_Coef(1:2,:); 0 0 0];
                        idx = idx+1;

                    case 'Specific Coordinate' %% Step 2 - Option 2: Convert to a particular Coor.

                        Specific_subMatrix = ((subMatrix)*R );
                        Specific_subDisp =  subDisp * R;
                        Specific_Corr_Coef = [Specific_subMatrix(:,1:2),ones(square_size^2, 1)]\Specific_subDisp;
                        % Save data
                        %coefficients( idx, :) = Local_Corr_Coef(1:3)';
                        coefficients{idx,1} = [Specific_Corr_Coef(1:2,:); 0 0 0];
                        idx = idx+1;


                end
            end
        end % end of UVW loop

    end
else
    coefficients = cell(N*M, 1);
    valid_idx = zeros(N*M, 1);

    parfor idx = 1:N*M
        [j, i] = ind2sub([M, N], idx);

        if task == 1 || task == 2 || task == 3
            %% task 1- 2 3 ------------
            R = [0.989437165332156,-0.014142182990485,0.144364060924945;-0.018226621689098,-0.999241912622796,0.033854241154337;0.143811981841816,-0.036271183027293,-0.988944849725807];
        elseif task == 4

            %%------------------
            %% task 4-------------
            R = [0.990006874274756	-0.00105708521970580	0.140817155786656;
-0.0122379091647181	-0.999442475470632	0.0328889703736955;
0.140487090040345	-0.0333709573727341	-0.989467857671003];
            %%------------------

        end


        if ~((i<1+Rad) || (i>N-Rad) || (j<1+Rad) || (j>M-Rad))
            %% Step 1 Get the dudx in Camera 1's Coor.
            % U/V/W = ax + by + cz +d ; Get a,b,c,d by Least squares method
            % Get the index of points in VSG
            Index = getSquareIDs(M, N, i, j, Rad, indexMap);
            % \[subDisp\] = \[subMatrix\] .\* \[Coef\]
            subMatrix = [Coordinate_ref(Index,:)]; % Lagrange strain: reference_coor
            subDisp = [Displacement{1,1},Displacement{1,2},Displacement{1,3}]; % Dimension: (2\*R+1)^2,1
            subDisp = subDisp(Index,:);

            switch Strain_Coordinate
                case 'Camera0'
                    % LSM
                    Camera0_Corr_Coef = [subMatrix,ones(square_size^2, 1)]\subDisp;
                    % Save data
                    coefficients{idx, 1} = Camera0_Corr_Coef(1:3);
                    valid_idx(idx) = 1;
                case 'Local' %% Step 2 - Option 1: Convert dudx to Local Coor.
                    % What was previously fitted was the four-dimensional
                    % surface of UVW, but now what needs to be re-projected
                    % is the three-dimensional surface of the actual shape. Z = ax + b y + c
                    subXY = [Coordinate_ref(Index,1:2),ones(square_size^2, 1)];
                    subZ = Coordinate_ref(Index,3);
                    Local_Corr_Coef = subXY\subZ;
                    % Rotation matrix
                    direction_z = [Local_Corr_Coef(1); Local_Corr_Coef(2) ;-1]; direction_z = direction_z/norm(direction_z);
                    direction_x = [1;0;0] - [1;0;0]'*direction_z*direction_z; direction_x = direction_x/norm(direction_x);
                    direction_y = cross(direction_z,direction_x);
                    R = [direction_x,direction_y,direction_z];
                    Local_subMatrix = ( (subMatrix ) * R);
                    Local_subDisp = subDisp * R;
                    Local_Corr_Coef = [Local_subMatrix(:,1:2),ones(square_size^2, 1)]\Local_subDisp;
                    % Save data
                    coefficients{idx,1} = [Local_Corr_Coef(1:2,:); 0 0 0];
                    valid_idx(idx) = 1;
                case 'Specific Coordinate' %% Step 2 - Option 2: Convert to a particular Coor.
                    Specific_subMatrix = ((subMatrix)*R );
                    Specific_subDisp = subDisp * R;
                    Specific_Corr_Coef = [Specific_subMatrix(:,1:2),ones(square_size^2, 1)]\Specific_subDisp;
                    % Save data
                    coefficients{idx,1} = [Specific_Corr_Coef(1:2,:); 0 0 0];
                    valid_idx(idx) = 1;
            end
        end
    end

    coefficients = coefficients(valid_idx == 1);

end
end

