function [] = plot_precondEffectS( generate_plots, ...
                                   len_xmesh, len_tmesh, n_sobolev_preconditioning_s_values, ...
                                   tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                                   initial_data_parameter_s, initial_data_parameter_a, ...
                                   sobolev_preconditioning_a ...
                                 )
  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.)
  optimization_parameter_defaults;

  % Set more default values for parameters
  if ~exist('generate_plots', 'var')
    generate_plots = true;
  end
  if ~exist('n_regularization_s_values', 'var')
    n_sobolev_preconditioning_s_values =50;
  end
   if ~exist('n_regularization_a_values', 'var')
    n_sobolev_preconditioning_a_values = 50;
   end

 

 
  output_headers = 'sobolev_preconditioning_a,J_final,||s_{final}-s_{true}||,||a_{final}-a_{true}||';
 

 
    % Choose range of sobolev_preconditioning_s values
    sobolev_preconditioning_a_values = linspace(1e-8, 10, n_sobolev_preconditioning_a_values)';

    % Initialize storage for output data
    J_values_final = zeros(size(sobolev_preconditioning_a_values))*NaN;
    dS_values_final = zeros(size(sobolev_preconditioning_a_values))*NaN;
    dA_values_final = zeros(size(sobolev_preconditioning_a_values))*NaN;

    % Discretization for both for forward and adjoint problem
    tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

    % Calculate corresponding s_true and a_true
    [~, ~, ~, ~, ~, s_true, ~, a_true] = true_solution(tmesh);

    % Since the mesh is fixed, we can compute the analytic solution once.
    s_true_values = s_true(tmesh)';
    a_true_values = a_true(tmesh)';

   
 
    for i = 1:n_sobolev_preconditioning_a_values
      sobolev_preconditioning_a = sobolev_preconditioning_a_values(i);
      try
        [jvals, svals, avals] = ...
          optimization( ...
                        len_xmesh, len_tmesh, ...
                        tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                        initial_data_parameter_s, initial_data_parameter_a, ...
                        sobolev_preconditioning_s, sobolev_preconditioning_a, ...
                        reconstruct_s, reconstruct_a, ...
                        false);

      J_values_final(i) = jvals(end);
      dS_values_final(i) = sqrt(trapz(tmesh, (svals(end, :) - s_true_values).^2));
      dA_values_final(i) = sqrt(trapz(tmesh, (avals(end, :) - a_true_values).^2));

      catch ME % Fail gracefully by ignoring the error
        warning('Error in optimization at i=%d (sobolev_preconditioning_s=%0.16f):', i, sobolev_preconditioning_s);
        disp(ME)
      end
    end % Loop over sobolev_preconditioning_s_values

    A(:,1)=sobolev_preconditioning_a_values;
    A(:,2)= J_values_final;
    A(:,3)= dS_values_final;
    A(:,4)= dA_values_final; 

    
    
    fid = fopen('Data_for_precond.txt','wt');
    fprintf(fid, ' Precond    J_values    S_Norm    A_Norm\n\n' );
    fclose(fid);
    
    dlmwrite('Data_for_precond.txt',A,'-append','delimiter','\t','newline', 'pc')

      %% Create plots
      % J(v_{final})
      subplot(3,1,1)
      plot(sobolev_preconditioning_a_values, J_values_final);
      xlabel('l_a')
      ylabel('J(v_{final})')

      % ||s_{final} - s_{true}||
      subplot(3,1,2)
      plot(sobolev_preconditioning_a_values, dS_values_final);
      xlabel('l_a')
      ylabel('||a_{final}-a_{true}||_{L_2}')

      % ||a_{final} - a_{true}||
      subplot(3,1,3)
      plot(sobolev_preconditioning_a_values, dA_values_final);
      xlabel('l_a')
      ylabel('||a_{final}-a_{true}||_{L_2}')
  end % if generate_plots

