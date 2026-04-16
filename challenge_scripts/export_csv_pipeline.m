function export_csv_pipeline(mat_file, frame_map, output_root)
%% Export 13-column submission CSVs from a coarse_test _strain.mat file.
%
% ALL frames go through the SAME compute path (disp interpolation + strain
% recomputation per VSG). No zero-shortcuts — noise frames and reference
% frame alike are computed, because noise strain IS a required output of
% the challenge (noise baseline evaluation).
%
% INPUTS:
%   mat_file    - Path to coarse_test_ws*_ss*_strain.mat
%   frame_map   - Struct array with fields:
%                   .mat_index    (integer, 1-based frame index in FinalResult)
%                   .csv_name     (char vector, e.g. '0200' -> frame_0200.csv)
%   output_root - Root folder containing VSG_23/, VSG_43/, VSG_63/ subdirs
%
% Produces:   output_root/VSG_XX/frame_YYYY.csv  for each (frame x VSG)
%
% VSG table (fixed): VSG 23 -> strain_size=2 (labeled 23, actual 26)
%                    VSG 43 -> strain_size=19 (exact)
%                    VSG 63 -> strain_size=39 (exact)

    fprintf('\n======= CSV EXPORT PIPELINE =======\n');
    fprintf('Loading %s ...\n', mat_file);
    S = load(mat_file);
    FinalResult = S.FinalResult;
    DICpara     = S.DICpara;
    RD_L        = S.RD_L;

    coords2D_sparse = RD_L.ResultFEMeshEachFrame{1,1}.coordinatesFEM;  % N x 2 [row,col]
    X_ref = FinalResult.CoordinatesNew{1,1};
    Y_ref = FinalResult.CoordinatesNew{1,2};
    Z_ref = FinalResult.CoordinatesNew{1,3};
    coords3D_ref_sparse = [X_ref, Y_ref, Z_ref];

    vsg_configs = {
        'VSG_23',  2;   % label 23, actual 26
        'VSG_43', 19;   % exact 43
        'VSG_63', 39;   % exact 63
    };
    num_vsg = size(vsg_configs, 1);
    num_frames = length(frame_map);

    % === Interpolate reference 3D coords ONCE (shared across all frames) ===
    zero_disp = zeros(size(coords3D_ref_sparse));
    [coords2D_dense, coords3D_dense, ~, mask_dense] = ...
        interpolate_to_dense_grid_from_mask(coords2D_sparse, coords3D_ref_sparse, ...
                                             zero_disp, DICpara, 1);
    fprintf('Dense grid: %d points, %d valid\n\n', ...
        size(coords2D_dense, 1), sum(mask_dense));

    pipeline_start = tic;

    % === Loop frames — SAME compute path for every frame ===
    for i = 1:num_frames
        fm = frame_map(i);
        fprintf('========== Frame %s [%d/%d] ==========\n', ...
                fm.csv_name, i, num_frames);

        % Extract sparse displacement for this frame (specimen coords)
        U = FinalResult.DisplacementNew{fm.mat_index, 1};
        V = FinalResult.DisplacementNew{fm.mat_index, 2};
        W = FinalResult.DisplacementNew{fm.mat_index, 3};
        disp_sparse_f = [U, V, W];

        u_mag = sqrt(U.^2 + V.^2 + W.^2);
        fprintf('  Sparse disp |U| mean=%.5f mm, max=%.5f mm\n', ...
            mean(u_mag, 'omitnan'), max(u_mag, [], 'omitnan'));

        % Interpolate displacement to dense grid
        [~, ~, disp_dense, ~] = interpolate_to_dense_grid_from_mask( ...
            coords2D_sparse, coords3D_ref_sparse, disp_sparse_f, DICpara, 1);

        % Recompute strain per VSG on the dense grid
        for v = 1:num_vsg
            vsg_name    = vsg_configs{v,1};
            strain_size = vsg_configs{v,2};

            fprintf('  --- Frame %s: %s [%d/%d] ---\n', ...
                fm.csv_name, vsg_name, v, num_vsg);

            % R = eye(3): CoordinatesNew is already in specimen coords (do NOT re-rotate)
            [exx, eyy, exy, e1, e2] = compute_strain_dense_grid_lowmem( ...
                coords2D_dense, coords3D_dense, disp_dense, ...
                strain_size, DICpara.winsize, 1, mask_dense, eye(3), 5000);

            out_dir = fullfile(output_root, vsg_name);
            export_challenge_csv(out_dir, fm.csv_name, ...
                coords2D_dense, coords3D_dense, disp_dense, ...
                exx, eyy, exy, e1, e2, mask_dense);
        end
        fprintf('  Frame %s done. Total elapsed: %.1f min\n\n', ...
            fm.csv_name, toc(pipeline_start)/60);
    end

    fprintf('======= CSV EXPORT COMPLETE (%.1f min) =======\n', ...
        toc(pipeline_start)/60);
end
