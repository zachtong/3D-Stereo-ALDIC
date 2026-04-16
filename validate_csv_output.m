%% Validate CSV export for Stereo DIC Challenge 2.1 Bespoke submission
% Standalone script to double-check CSV outputs produced by run_export_csv.m
% are structurally correct and physically plausible BEFORE submission.
%
% Read-only: does not modify any files.
% Run from the project root. Requires results/csv_output/VSG_{23,43,63}/ to exist.

close all; clear; clc;

all_passed = true;
warnings_count = 0;

expected_cols = {'X_image','Y_image','X','Y','Z','U','V','W', ...
                 'epsilon_xx','epsilon_yy','epsilon_xy','epsilon_1','epsilon_2'};

base_dir = fullfile('results', 'csv_output');

%% Check 1: CSV count per VSG directory
fprintf('\n=== Check 1: CSV count per VSG directory ===\n');
try
    vsg_list = [23, 43, 63];
    for k = 1:numel(vsg_list)
        vsg = vsg_list(k);
        dir_path = fullfile(base_dir, sprintf('VSG_%d', vsg));
        if exist(dir_path, 'dir') ~= 7
            error('Directory missing: %s', dir_path);
        end
        files = dir(fullfile(dir_path, 'frame_*.csv'));
        n = numel(files);
        fprintf('VSG_%d: %d CSVs\n', vsg, n);
        assert(n == 18, sprintf('VSG_%d: expected 18 CSV files, found %d', vsg, n));
    end
    fprintf('Check 1 PASSED\n');
catch ME
    fprintf(2, 'Check 1 FAILED: %s\n', ME.message);
    all_passed = false;
end

%% Check 2: CSV structure (schema)
fprintf('\n=== Check 2: CSV structure (schema) ===\n');
try
    csv_file = fullfile(base_dir, 'VSG_43', 'frame_0200.csv');
    if exist(csv_file, 'file') ~= 2
        error('CSV file missing: %s', csv_file);
    end
    T = readtable(csv_file);
    assert(width(T) == 13, sprintf('Expected 13 columns, got %d', width(T)));
    names = T.Properties.VariableNames;
    assert(numel(names) == numel(expected_cols) && ...
           all(strcmp(names, expected_cols)), ...
           'Column names/order mismatch');
    fprintf('frame_0200 rows=%d, disp U range=[%.3f, %.3f] mm\n', ...
            height(T), min(T.U), max(T.U));
    fprintf('Check 2 PASSED\n');
catch ME
    fprintf(2, 'Check 2 FAILED: %s\n', ME.message);
    all_passed = false;
end

%% Check 3: Reference frame (frame_0000) is effectively zero
fprintf('\n=== Check 3: Reference frame (frame_0000) effectively zero ===\n');
try
    csv_file = fullfile(base_dir, 'VSG_43', 'frame_0000.csv');
    if exist(csv_file, 'file') ~= 2
        error('CSV file missing: %s', csv_file);
    end
    T_ref = readtable(csv_file);
    u_mag = sqrt(T_ref.U.^2 + T_ref.V.^2 + T_ref.W.^2);
    u_mag_max = max(u_mag, [], 'omitnan');
    exx_max = max(abs(T_ref.epsilon_xx), [], 'omitnan');
    fprintf('frame_0000 (reference): |U| max=%.6f mm, |exx| max=%.6f\n', ...
            u_mag_max, exx_max);
    assert(u_mag_max < 1e-6, 'Reference disp should be ~0 (self-correlation)');
    assert(exx_max < 1e-6, 'Reference strain should be ~0');
    fprintf('\n*** frame_0000 is computed and effectively zero ***\n');
    fprintf('Check 3 PASSED\n');
catch ME
    fprintf(2, 'Check 3 FAILED: %s\n', ME.message);
    all_passed = false;
end

%% Check 4: Noise frames have small displacement
fprintf('\n=== Check 4: Noise frame has small displacement ===\n');
try
    csv_file = fullfile(base_dir, 'VSG_43', 'frame_0003.csv');
    if exist(csv_file, 'file') ~= 2
        error('CSV file missing: %s', csv_file);
    end
    T_noise = readtable(csv_file);
    u_mag = sqrt(T_noise.U.^2 + T_noise.V.^2 + T_noise.W.^2);
    u_mean = mean(u_mag, 'omitnan');
    u_max = max(u_mag, [], 'omitnan');
    fprintf('Noise frame_0003: |U| mean=%.4f mm, max=%.4f mm\n', u_mean, u_max);
    assert(u_mean < 0.02, ...
           'Noise displacement should be << loading displacement (> 0.02 mm means something is wrong)');
    fprintf('Check 4 PASSED\n');
catch ME
    fprintf(2, 'Check 4 FAILED: %s\n', ME.message);
    all_passed = false;
end

%% Check 5: Loading frames have monotonically growing displacement
fprintf('\n=== Check 5: Loading frames monotonicity ===\n');
try
    csvs = {'0200','0600','1000','1400','1800','2200','2600'};
    u_means = zeros(numel(csvs), 1);
    for i = 1:numel(csvs)
        csv_file = fullfile(base_dir, 'VSG_43', sprintf('frame_%s.csv', csvs{i}));
        if exist(csv_file, 'file') ~= 2
            error('CSV file missing: %s', csv_file);
        end
        T = readtable(csv_file);
        u_means(i) = mean(sqrt(T.U.^2 + T.V.^2 + T.W.^2), 'omitnan');
    end
    disp(table(csvs', u_means, 'VariableNames', {'frame', 'mean_disp_mm'}));

    % Lenient monotonicity check: flag only significant decreases (>0.01 mm)
    violations = {};
    for i = 1:numel(csvs)-1
        if u_means(i+1) < u_means(i) - 0.01
            violations{end+1} = sprintf(...
                'frame %s -> %s: %.4f -> %.4f (drop of %.4f mm)', ...
                csvs{i}, csvs{i+1}, u_means(i), u_means(i+1), ...
                u_means(i) - u_means(i+1)); %#ok<AGROW>
        end
    end

    if isempty(violations)
        fprintf('Monotonicity: OK\n');
        fprintf('Check 5 PASSED\n');
    else
        fprintf(2, 'Monotonicity WARNING: %d violation(s) found:\n', numel(violations));
        for v = 1:numel(violations)
            fprintf(2, '  - %s\n', violations{v});
        end
        fprintf(2, '(Monotonicity is an approximate check; not fatal.)\n');
        warnings_count = warnings_count + numel(violations);
        fprintf('Check 5 PASSED (with %d monotonicity warning(s))\n', numel(violations));
    end
catch ME
    fprintf(2, 'Check 5 FAILED: %s\n', ME.message);
    all_passed = false;
end

%% Final summary
if all_passed
    if warnings_count == 0
        fprintf('\n✓✓✓ ALL CHECKS PASSED — CSV submission is ready ✓✓✓\n');
    else
        fprintf('\n✓✓✓ ALL CHECKS PASSED with %d non-fatal warning(s) — review above ✓✓✓\n', ...
                warnings_count);
    end
else
    error('validate_csv_output:failures', ...
          'Some checks failed; see messages above.');
end
