function [coordinates2D_dense, coordinates3D_dense, displacement_dense, mask_dense] = ...
    interpolate_to_dense_grid_from_mask(coordinates2D_original, coordinates3D_original, ...
    displacement_original, DICpara, target_stepsize)
%% Interpolate displacement field to a dense grid using ImgRefMask
% This version extracts the ROI directly from DICpara.ImgRefMask to ensure
% consistent grid across all frames
%
% INPUT:
%   coordinates2D_original - Original 2D coordinates (N x 2) [row, col] in pixels
%   coordinates3D_original - Original 3D coordinates (N x 3) [X, Y, Z] in mm
%   displacement_original - Original displacement (N x 3) [U, V, W] in mm
%   DICpara - DIC parameters structure (must contain ImgRefMask)
%   target_stepsize - Target step size (typically 1 pixel)
%
% OUTPUT:
%   coordinates2D_dense - Dense 2D grid coordinates (M x 2) [row, col]
%   coordinates3D_dense - Interpolated 3D coordinates
%   displacement_dense - Interpolated displacement field
%   mask_dense - Valid points mask (1 = valid, 0 = outside ROI)

% Extract ROI from ImgRefMask
if ~isfield(DICpara, 'ImgRefMask')
    error('DICpara.ImgRefMask not found! Cannot determine ROI.');
end

mask = DICpara.ImgRefMask;

% Find ROI bounding box from mask
[mask_rows, mask_cols] = find(mask);
if isempty(mask_rows)
    error('ImgRefMask is empty! Cannot determine ROI.');
end

% ROI bounds in [row, col] convention
row_min = min(mask_rows);
row_max = max(mask_rows);
col_min = min(mask_cols);
col_max = max(mask_cols);

% Extract original 2D coordinates (coordinatesFEM format: [row, col])
row_orig = coordinates2D_original(:, 1);
col_orig = coordinates2D_original(:, 2);

% Create dense grid with target_stepsize using ROI bounds from mask
% meshgrid(col_range, row_range) produces:
%   col_dense varies along columns, row_dense varies along rows
[col_dense, row_dense] = meshgrid(col_min:target_stepsize:col_max, ...
                                   row_min:target_stepsize:row_max);

% Flatten the dense grid for interpolation
row_dense_vec = row_dense(:);
col_dense_vec = col_dense(:);
coordinates2D_dense = [row_dense_vec, col_dense_vec];

%% Interpolate 3D coordinates
% Create scattered interpolant for each coordinate
% Training order: (row, col) - must match query order below
F_X = scatteredInterpolant(row_orig, col_orig, coordinates3D_original(:, 1), 'linear', 'none');
F_Y = scatteredInterpolant(row_orig, col_orig, coordinates3D_original(:, 2), 'linear', 'none');
F_Z = scatteredInterpolant(row_orig, col_orig, coordinates3D_original(:, 3), 'linear', 'none');

% Interpolate - query with SAME (row, col) order as training
X_dense = F_X(row_dense_vec, col_dense_vec);
Y_dense = F_Y(row_dense_vec, col_dense_vec);
Z_dense = F_Z(row_dense_vec, col_dense_vec);

coordinates3D_dense = [X_dense, Y_dense, Z_dense];

%% Interpolate displacement field
% Create scattered interpolant for each displacement component
% Training order: (row, col) - must match query order below
F_U = scatteredInterpolant(row_orig, col_orig, displacement_original(:, 1), 'linear', 'none');
F_V = scatteredInterpolant(row_orig, col_orig, displacement_original(:, 2), 'linear', 'none');
F_W = scatteredInterpolant(row_orig, col_orig, displacement_original(:, 3), 'linear', 'none');

% Interpolate - query with SAME (row, col) order as training
U_dense = F_U(row_dense_vec, col_dense_vec);
V_dense = F_V(row_dense_vec, col_dense_vec);
W_dense = F_W(row_dense_vec, col_dense_vec);

displacement_dense = [U_dense, V_dense, W_dense];

%% Create mask for valid points
% Method 1: Points with NaN values from interpolation are outside data coverage
interp_mask = ~isnan(X_dense) & ~isnan(Y_dense) & ~isnan(Z_dense) & ...
              ~isnan(U_dense) & ~isnan(V_dense) & ~isnan(W_dense);

% Method 2: Check against ImgRefMask
% Convert to subscripts and check mask
mask_from_ref = false(size(row_dense_vec));
for i = 1:length(row_dense_vec)
    ri = round(row_dense_vec(i));
    ci = round(col_dense_vec(i));

    % Check if within image bounds: mask is (num_rows x num_cols)
    if ri >= 1 && ri <= size(mask, 1) && ci >= 1 && ci <= size(mask, 2)
        mask_from_ref(i) = mask(ri, ci) > 0;
    end
end

% Combined mask: both interpolation successful AND within original mask
mask_dense = interp_mask & mask_from_ref;

end
