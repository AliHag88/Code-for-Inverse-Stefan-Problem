function [] = plot_precondEffectS( generate_data, generate_plots, ...
  len_xmesh, len_tmesh, n_regularization_s_values, ...
  tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
  initial_data_parameter_s, initial_data_parameter_a, ...
  regularization_a ...
  )
  % Set default values of input parameters
  if ~exist('generate_data', 'var')
    generate_data = true;
  end
  if ~exist('generate_plots', 'var')
    generate_plots = false;
  end
  if ~exist('len_xmesh', 'var')
    len_xmesh = 100;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 40;
  end
  if ~exist('tolerance', 'var')
    tolerance = 1e-9;
  end
  
  % Flags below change the optimization routine
  if ~exist('num_iterations', 'var')
    num_iterations = 40;
  end
  if ~exist('num_sub_iterations', 'var')
    num_sub_iterations = 20;
  end
  
  % See initial_setup.m for the following parameters
  if ~exist('use_synthetic_data', 'var')
    use_synthetic_data = true;
  end
  if ~exist('initial_data_parameter_s', 'var')
    initial_data_parameter_s = 0.6;
  end
  if ~exist('initial_data_parameter_a', 'var')
    initial_data_parameter_a = 0.6;
  end
  % Preconditioning parameters 
  if ~exist('regularization_a', 'var')
    regularization_a = 0.22;
  end
  
  output_filename = sprintf('precondEffectS_Nx%d_Nt%d.csv', len_xmesh, len_tmesh);
  output_headers = 'regularization_s,J_final,||s_{final}-s_{true}||,||a_{final}-a_{true}||';
  output_format = '%.16f,';
  output_size = [1, 3];
  
  % Generate example data
  if generate_data
    regularization_s_values = linspace(1e-8,1,n_regularization_s_values)';
    
    J_values_final = zeros(size(regularization_s_values))*NaN;
    dS_values_final = zeros(size(regularization_s_values))*NaN;
    dA_values_final = zeros(size(regularization_s_values))*NaN;
    
    t_final = 1;
    
    % Discretization for both for forward and adjoint problem
    tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)
    
    [~, ~, ~, ~, ~, s_true, ~, a_true] = true_solution(tmesh);
    
    s_true_values = s_true(tmesh)';
    a_true_values = a_true(tmesh)';
    
    for i = 1:n_regularization_s_values
      regularization_s = regularization_s_values(i);
      [jvals, svals, avals] = optimization( ...
          len_xmesh, len_tmesh, ...
          tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
          initial_data_parameter_s, initial_data_parameter_a, ...
          regularization_s, regularization_a, false);
      
      J_values_final(i) = jvals(end);
      dS_values_final(i) = sqrt(trapz(tmesh, (svals(end, :) - s_true_values).^2));
      dA_values_final(i) = sqrt(trapz(tmesh, (avals(end, :) - a_true_values).^2));      
    end % Loop over regularization_s_values
    
    % Output data to file in CSV format
    fileID = fopen(output_filename, 'w');
    
    fprintf(fileID, output_headers);
    for i = 1:length(n_regularization_s_values)
      fprintf(fileID, [repmat(output_format, output_size) '\n'], ...
          regularization_s_values(i), ...
          J_values_final(i), ...
          dS_values_final(i), ...
          dA_values_final(i) ...
          );
    end % For loop
    fclose(fileID);
  end % if generate_data
  
  if generate_plots
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
            if tline == -1 % Documented as EOL marker
                break
            end
            values_tmp = sscanf(tline, '%f,', output_size);
            regularization_s_values(end + 1) = values_tmp(1);
            J_values_final(end + 1) = values_tmp(2);
            dS_values_final(end + 1) = values_tmp(3);
            dA_values_final(end + 1) = values_tmp(4);
          end % while loop
      end % if variables don't exist
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
      
      subplot(3,1,3)
      plot(regularization_s_values, dA_values_final);
      xlabel('l_s')
      ylabel('||a_{final}-a_{true}||_{L_2}')
  end % if generate_plots
end
