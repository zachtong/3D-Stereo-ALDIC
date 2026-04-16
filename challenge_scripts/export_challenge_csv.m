function export_challenge_csv(output_folder, frame_name, ...
    coordinates2D, coordinates3D, displacement, ...
    strain_exx, strain_eyy, strain_exy, strain_e1, strain_e2, mask)
%% Export results to CSV format for Stereo DIC Challenge 2.1
% This function exports the results in the required format:
% 13 columns: X_image, Y_image, X, Y, Z, U, V, W, εxx, εyy, εxy, ε1, ε2
%
% INPUT:
%   output_folder - Output folder path
%   frame_name - Frame name (e.g., '0000', '0001', ...)
%   coordinates2D - 2D pixel coordinates (N × 2)
%   coordinates3D - 3D world coordinates (N × 3) [X, Y, Z] in mm
%   displacement - Displacement field (N × 3) [U, V, W] in mm
%   strain_exx, strain_eyy, strain_exy - Strain components (N × 1)
%   strain_e1, strain_e2 - Principal strains (N × 1)
%   mask - Valid points mask (N × 1), only export valid points

% Create output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Filter valid points
valid_idx = find(mask & ~isnan(strain_exx));
num_valid = length(valid_idx);

if num_valid == 0
    warning('No valid points to export for frame %s', frame_name);
    return;
end

% Extract valid data
% coordinates2D format is [row, col] from coordinatesFEM
% Challenge convention: X_image = x_px = row (0-indexed), Y_image = y_px = col (0-indexed)
X_image = coordinates2D(valid_idx, 1) - 1;  % row -> X_image (0-indexed)
Y_image = coordinates2D(valid_idx, 2) - 1;  % col -> Y_image (0-indexed)
X = coordinates3D(valid_idx, 1);
Y = coordinates3D(valid_idx, 2);
Z = coordinates3D(valid_idx, 3);
U = displacement(valid_idx, 1);
V = displacement(valid_idx, 2);
W = displacement(valid_idx, 3);
exx = strain_exx(valid_idx);
eyy = strain_eyy(valid_idx);
exy = strain_exy(valid_idx);
e1 = strain_e1(valid_idx);
e2 = strain_e2(valid_idx);

% Create output table
output_table = table(X_image, Y_image, X, Y, Z, U, V, W, exx, eyy, exy, e1, e2, ...
    'VariableNames', {'X_image', 'Y_image', 'X', 'Y', 'Z', 'U', 'V', 'W', ...
                      'epsilon_xx', 'epsilon_yy', 'epsilon_xy', 'epsilon_1', 'epsilon_2'});

% Output file path
output_file = fullfile(output_folder, ['frame_', frame_name, '.csv']);

% Write to CSV
writetable(output_table, output_file);

fprintf('      wrote %s (%d pts)\n', output_file, num_valid);

end
