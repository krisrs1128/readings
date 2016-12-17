function [] = epca_wrapper(input_data, output_dir, path_to_epca, exp_fam_type, r, sigma_sq, bino_n)

  % Wrapper for Liu & Dobriban Exponential Family PCA
  %
  % This calls exp_fam_pca on data stored in a csv file and then saves the
  % results to a user specified path. The original function is defined here
  % https://github.com/lydiatliu/epca/blob/master/software/exp_fam_pca.m
  %
  % Input:
  %   input_data - The full path and filename for the csv file containing
  %     the data on which to apply epca.
  %   output_dir - The directory to save the resulting pca data.
  %   path_to_epca - The full path to the epca/software/ directory
  % Output:
  %   None.
  % Side Effects:
  %   Writes files V.txt, U.txt, and eval.txt to the output_dir directory.
  Y = csvread(input_data);
  addpath([path_to_epca]);
  [~, ~, white_V, white_U, white_eval, ~, ~, ~, ~] = exp_fam_pca(Y, exp_fam_type, r, sigma_sq, bino_n);

  dlmwrite(strcat(output_dir, '/V.txt'), white_V);
  dlmwrite(strcat(output_dir, '/U.txt'), white_U);
  dlmwrite(strcat(output_dir, '/evals.txt'), white_eval);
