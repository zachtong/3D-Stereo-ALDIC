%% Driver: Export all 54 challenge CSVs from coarse test + noise results
close all; clear; clc;

addpath('./challenge_scripts');
addpath('./func');
addpath('./func_quadtree');

output_root = 'results/csv_output';

%% --- Frame map for LOADING + FAILURE dataset ---
% Loading mat file has 136 frames; index = XXXX/20 + 1
loading_map = struct('mat_index', {}, 'csv_name', {});
loading_map(end+1) = struct('mat_index', 1,   'csv_name', '0000');  % reference (self-corr)
loading_map(end+1) = struct('mat_index', 11,  'csv_name', '0200');
loading_map(end+1) = struct('mat_index', 31,  'csv_name', '0600');
loading_map(end+1) = struct('mat_index', 51,  'csv_name', '1000');
loading_map(end+1) = struct('mat_index', 71,  'csv_name', '1400');
loading_map(end+1) = struct('mat_index', 91,  'csv_name', '1800');
loading_map(end+1) = struct('mat_index', 111, 'csv_name', '2200');
loading_map(end+1) = struct('mat_index', 131, 'csv_name', '2600');
loading_map(end+1) = struct('mat_index', 132, 'csv_name', '2620');
loading_map(end+1) = struct('mat_index', 133, 'csv_name', '2640');
loading_map(end+1) = struct('mat_index', 134, 'csv_name', '2660');
loading_map(end+1) = struct('mat_index', 135, 'csv_name', '2680');
loading_map(end+1) = struct('mat_index', 136, 'csv_name', '2700');

%% --- Frame map for NOISE dataset ---
% Noise mat file has 6 frames; f=2..6 -> CSV 0001..0005
% (f=1 in noise dataset is also reference 0000 — we use Loading's instead
%  to have ONE canonical reference in submission)
noise_map = struct('mat_index', {}, 'csv_name', {});
noise_map(end+1) = struct('mat_index', 2, 'csv_name', '0001');
noise_map(end+1) = struct('mat_index', 3, 'csv_name', '0002');
noise_map(end+1) = struct('mat_index', 4, 'csv_name', '0003');
noise_map(end+1) = struct('mat_index', 5, 'csv_name', '0004');
noise_map(end+1) = struct('mat_index', 6, 'csv_name', '0005');

%% --- Execute ---
loading_mat = 'results/loading_and_failure/coarse_test_ws24_ss8_strain.mat';
noise_mat   = 'results/noise/coarse_test_noise_ws24_ss8_strain.mat';

fprintf('\n>>>>> LOADING + FAILURE DATASET <<<<<\n');
export_csv_pipeline(loading_mat, loading_map, output_root);

fprintf('\n>>>>> NOISE DATASET <<<<<\n');
export_csv_pipeline(noise_mat, noise_map, output_root);

fprintf('\nALL DONE. Expect 54 CSVs in %s/VSG_{23,43,63}/\n', output_root);
