function [] = plot_precondEffectS( generate_data, generate_plots, ...
                                   len_xmesh, len_tmesh, n_regularization_s_values, ...
                                   tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                                   initial_data_parameter_s, initial_data_parameter_a, ...
                                   regularization_a ...
                                 )
  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.
  optimization_parameter_defaults;

  % Set more default values for parameters
  if ~exist('generate_data', 'var')
    generate_data = true;
  end
  if ~exist('generate_plots', 'var')
    generate_plots = false;
  end
  if ~exist('n_regularization_s_values', 'var')
    n_regularization_s_values = 10;
  end
  
  % Compute filename corresponding to this example
  output_filename = sprintf('precondEffectS_Nx%d_Nt%d.csv', len_xmesh, len_tmesh);
  % Headers and file format for output
  output_headers = 'regularization_s,J_final,||s_{final}-s_{true}||,||a_{final}-a_{true}||';
  output_format = '%.16f,';
  output_size = [1, 3];
  
  % Generate example data
  if generate_data
    % Choose range of regularization_s values
    regularization_s_values = linspace(1e-8,1,n_regularization_s_values)';

    % Initialize storage for output data
    J_values_final = zeros(size(regularization_s_values))*NaN;
    dS_values_final = zeros(size(regularization_s_values))*NaN;
    dA_values_final = zeros(size(regularization_s_values))*NaN;

    % Set final moment (as in optimization.m code)
    t_final = 1;

    % Discretization for both for forward and adjoint problem
    tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

    % Calculate corresponding s_true and a_true
    [~, ~, ~, ~, ~, s_true, ~, a_true] = true_solution(tmesh);

    % Since the mesh is fixed, we can compute the analytic solution once.
    s_true_values = s_true(tmesh)';
    a_true_values = a_true(tmesh)';

    % Open file to export data in CSV format
    fileID = fopen(output_filename, 'w');

    % Output headers to file
    fprintf(fileID, output_headers);
    for i = 1:n_regularization_s_values
      regularization_s = regularization_s_values(i);
      try
        [jvals, svals, avals] = ...
          optimization( ...
                        len_xmesh, len_tmesh, ...
                        tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                        initial_data_parameter_s, initial_data_parameter_a, ...
                        regularization_s, regularization_a, false);

      J_values_final(i) = jvals(end);
      dS_values_final(i) = sqrt(trapz(tmesh, (svals(end, :) - s_true_values).^2));
      dA_values_final(i) = sqrt(trapz(tmesh, (avals(end, :) - a_true_values).^2));

      catch ME # Fail gracefully by ignoring the error
        warning(sprintf('Error in optimization at i=%d (regularization_s=%0.16f):', i, regularization_s));
        disp(ME)
      end

      % Output to file
      fprintf(fileID, [repmat(output_format, output_size) '\n'], ...
              regularization_s_values(i), ...
              J_values_final(i), ...
              dS_values_final(i), ...
              dA_values_final(i) ...
             );

    end % Loop over regularization_s_values

    % Close output file
    fclose(fileID);

  end % if generate_data
  
  if generate_plots
    % Check if the output values exist. If the don't, get them from the
    % output file.
      if ~exist('J_values_final', 'var') && isfile(output_filename)
        fileID = fopen(output_filename, 'r');
        fgetl(fileID); % Skip header

        % Initialize storage
        regularization_s_values = [];
        J_values_final = [];
        dS_values_final = [];
        dA_values_final = [];

        % Read all values from file
        while true
          tline = fgetl(fileID);
          if tline == -1 % Documented as EOL marker in fgetl
            break
          end
          % Parse output line
          values_tmp = sscanf(tline, '%f,', output_size);
          regularization_s_values(end + 1) = values_tmp(1);
          J_values_final(end + 1) = values_tmp(2);
          dS_values_final(end + 1) = values_tmp(3);
          dA_values_final(end + 1) = values_tmp(4);
        end % while loop
      end % if variables don't exist

      %% Create plots
      % J(v_{final})
      subplot(3,1,1)
      plot(regularization_s_values, J_values_final);
      xlabel('l_s')
      ylabel('J(v_{final})')

      % ||s_{final} - s_{true}||
      subplot(3,1,2)
      plot(regularization_s_values, dS_values_final);
      xlabel('l_s')
      ylabel('||s_{final}-s_{true}||_{L_2}')

      % ||a_{final} - a_{true}||
      subplot(3,1,3)
      plot(regularization_s_values, dA_values_final);
      xlabel('l_s')
      ylabel('||a_{final}-a_{true}||_{L_2}')
  end % if generate_plots
end