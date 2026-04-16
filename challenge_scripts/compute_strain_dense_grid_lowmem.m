function [strain_exx, strain_eyy, strain_exy, strain_e1, strain_e2] = ...
    compute_strain_dense_grid_lowmem(coordinates2D, coordinates3D, displacement, ...
    strain_size, winsize, stepsize, mask, coord_transform_R, batch_size)
%% Low-memory version: Compute strain on dense grid with batching
% This function processes points in batches to reduce memory usage
%
% Additional parameter:
%   batch_size - Number of points to process at once (default: 5000)

if nargin < 9
    batch_size = 5000;  % Default batch size
end

% Calculate VSG size
VSG_size = (strain_size - 1) * stepsize + 1 + winsize;

% Filter valid points
valid_idx = find(mask);
num_valid_points = length(valid_idx);

% Initialize strain arrays
strain_exx = nan(size(mask));
strain_eyy = nan(size(mask));
strain_exy = nan(size(mask));
strain_e1 = nan(size(mask));
strain_e2 = nan(size(mask));

% Calculate strain_length for neighbor search
strain_length = (strain_size - 1) * stepsize + 1;

% Get valid point data
valid_coords2D = coordinates2D(valid_idx, :);
valid_coords3D = coordinates3D(valid_idx, :);
valid_disp = displacement(valid_idx, :);

% Number of neighbors to search
Knn_number = min(2 * strain_size^2, num_valid_points);
if Knn_number < 9
    Knn_number = 9;
end

% Calculate number of batches
num_batches = ceil(num_valid_points / batch_size);
fprintf('    VSG=%d px, K=%d, %d valid pts, %d batches\n', ...
        VSG_size, Knn_number, num_valid_points, num_batches);

overall_start = tic;
next_progress_pct = 10;  % Print at 10%, 20%, ..., 100%

for batch = 1:num_batches
    % Determine batch range
    start_idx = (batch - 1) * batch_size + 1;
    end_idx = min(batch * batch_size, num_valid_points);
    batch_indices = start_idx:end_idx;

    % KNN search for this batch (silent)
    [neighborInd, distance] = knnsearch(valid_coords2D, valid_coords2D(batch_indices, :), ...
        'K', Knn_number, 'Distance', 'chebychev');

    % Process each point in this batch
    for local_i = 1:length(batch_indices)
        global_i = batch_indices(local_i);

        % Find points within VSG
        OutsideOrNot = distance(local_i, :) > 0.5 * strain_length;
        numInsideVSG = find(OutsideOrNot, 1, 'first') - 1;

        if isempty(numInsideVSG)
            numInsideVSG = Knn_number;
        end

        % Need at least 9 points
        if numInsideVSG < 9
            continue;
        end

        % Get neighbor indices
        tempIndex = neighborInd(local_i, 1:numInsideVSG);

        % Get 3D coordinates and displacements
        subMatrix = valid_coords3D(tempIndex, :);
        subDisp = valid_disp(tempIndex, :);

        % Compute local coordinate system if needed
        if ~isempty(coord_transform_R)
            R = coord_transform_R;
        else
            subXY = [subMatrix(:, 1:2), ones(numInsideVSG, 1)];
            subZ = subMatrix(:, 3);
            Local_Corr_Coef = subXY \ subZ;

            direction_z = [Local_Corr_Coef(1); Local_Corr_Coef(2); -1];
            direction_z = direction_z / norm(direction_z);
            direction_x = [1; 0; 0] - [1; 0; 0]' * direction_z * direction_z;
            direction_x = direction_x / norm(direction_x);
            direction_y = cross(direction_z, direction_x);
            R = [direction_x, direction_y, direction_z];
        end

        % Transform to local coordinate system
        Local_subMatrix = subMatrix * R;
        Local_subDisp = subDisp * R;

        % Fit displacement gradients
        try
            Local_Corr_Coef = [Local_subMatrix(:, 1:2), ones(numInsideVSG, 1)] \ Local_subDisp;

            % Displacement gradient components
            u_x = Local_Corr_Coef(1, 1);
            u_y = Local_Corr_Coef(2, 1);
            v_x = Local_Corr_Coef(1, 2);
            v_y = Local_Corr_Coef(2, 2);
            w_x = Local_Corr_Coef(1, 3);
            w_y = Local_Corr_Coef(2, 3);

            % Deformation gradient tensor
            F = [1+u_x, u_y, 0;
                 v_x, 1+v_y, 0;
                 w_x, w_y, 1];

            % Green-Lagrangian strain tensor: E = 0.5*(F'*F - I)
            E = 0.5 * (F' * F - eye(3));
            exx = E(1, 1);
            eyy = E(2, 2);
            exy = E(1, 2);

            % Principal strains
            strain_mean = (exx + eyy) / 2;
            strain_diff = (exx - eyy) / 2;
            strain_radius = sqrt(strain_diff^2 + exy^2);

            e1 = strain_mean + strain_radius;
            e2 = strain_mean - strain_radius;

            % Store results
            global_idx = valid_idx(global_i);
            strain_exx(global_idx) = exx;
            strain_eyy(global_idx) = eyy;
            strain_exy(global_idx) = exy;
            strain_e1(global_idx) = e1;
            strain_e2(global_idx) = e2;
        catch
            continue;
        end
    end

    % Progress report at 10% intervals
    pct = 100 * end_idx / num_valid_points;
    if pct >= next_progress_pct || batch == num_batches
        elapsed = toc(overall_start);
        rate = end_idx / elapsed;
        if end_idx < num_valid_points
            eta_str = sprintf('ETA %.1f min', (num_valid_points - end_idx) / rate / 60);
        else
            eta_str = 'done';
        end
        fprintf('    %3.0f%% | %.1f min | %s\n', pct, elapsed/60, eta_str);
        next_progress_pct = floor(pct/10)*10 + 10;
    end

    % Clear batch variables
    clear neighborInd distance;
end

fprintf('    Strain points: %d valid / %d total, %.1f min\n', ...
    sum(~isnan(strain_exx)), length(strain_exx), toc(overall_start)/60);

end
