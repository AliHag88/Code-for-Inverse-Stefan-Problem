function [] = plot_initialApproachAS( generate_data, generate_plots, ...
                                      len_xmesh, len_tmesh, ...
                                      n_initial_data_parameter_s_values, n_initial_data_parameter_a_values, ...
                                      tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                                      sobolev_preconditioning_s, sobolev_preconditioning_a ...
                                 )
  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.)
  optimization_parameter_defaults;

  % Set more default values for parameters
  if ~exist('generate_data', 'var')
    generate_data = true;
  end
  if ~exist('generate_plots', 'var')
    generate_plots = false;
  end
  if ~exist('n_initial_data_parameter_s_values', 'var')
    n_initial_data_parameter_s_values = 10;
  end
  if ~exist('n_initial_data_parameter_a_values', 'var')
    n_initial_data_parameter_a_values = 10;
  end

  % Compute filename corresponding to this example
  output_filename = sprintf('initialApproachAS_Nx%d_Nt%d.csv', len_xmesh, len_tmesh);

  plot_filename_2d = sprintf('initialApproachAS_2D_Nx%d_Nt%d.pdf', len_xmesh, len_tmesh);
  plot_filename_1d = sprintf('initialApproachAS_1D_Nx%d_Nt%d.pdf', len_xmesh, len_tmesh);

  % Headers and file format for output
  output_headers = 'initial_data_parameter_s,initial_data_parameter_a,J_final,||s_{final}-s_{true}||,||a_{final}-a_{true}||\n';
  output_format = '%.16f,';
  output_size = [1, 5];

  % Generate example data
  if generate_data
    % Choose range of initial_data_parameter_s values
    initial_data_parameter_s_values = linspace(-1, 1, n_initial_data_parameter_s_values)';
    % Choose range of initial_data_parameter_a values
    initial_data_parameter_a_values = linspace(-1, 1, n_initial_data_parameter_a_values);
    
    % Initialize storage for output data
    J_values_final = zeros(n_initial_data_parameter_s_values, n_initial_data_parameter_a_values)*NaN;
    dS_values_final = zeros(size(J_values_final))*NaN;
    dA_values_final = zeros(size(J_values_final))*NaN;
    
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
    % For each row
    for i = 1:n_initial_data_parameter_s_values
        initial_data_parameter_s = initial_data_parameter_s_values(i);
        % For each column
        for j = 1:n_initial_data_parameter_a_values
            initial_data_parameter_a = initial_data_parameter_a_values(j);

            try
                % Run optimization algorithm with chosen parameters
                [jvals, svals, avals] = ...
                    optimization( ...
                    len_xmesh, len_tmesh, ...
                    tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                    initial_data_parameter_s, initial_data_parameter_a, ...
                    sobolev_preconditioning_s, sobolev_preconditioning_a, ...
                    reconstruct_s, reconstruct_a, ...
                    false);

                J_values_final(i, j) = jvals(end);
                dS_values_final(i, j) = sqrt(trapz(tmesh, (svals(end, :) - s_true_values).^2));
                dA_values_final(i, j) = sqrt(trapz(tmesh, (avals(end, :) - a_true_values).^2));

            catch ME % Fail gracefully by ignoring the error
                warning('Error in optimization at (i,j)=(%d,%d) (initial_data_parameter_s=%0.16f, initial_data_parameter_a=%0.16f):', i, j, initial_data_parameter_s, initial_data_parameter_a);
                disp(ME)
            end

            % Output to file
            fprintf(fileID, [repmat(output_format, output_size) '\n'], ...
                initial_data_parameter_s, ...
                initial_data_parameter_a, ...
                J_values_final(i, j), ...
                dS_values_final(i, j), ...
                dA_values_final(i, j) ...
                );

        end % Loop over initial_data_parameter_a_values
        % Write another newline after each column
        fprintf(fileID, '\n');
    end % Loop over initial_data_parameter_s_values

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
        IDS = [];
        IDA = [];
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
            if isempty(values_tmp)
                warning('Failed to parse line "%s"', tline);
                continue
            end
            IDS(end + 1) = values_tmp(1);
            IDA(end + 1) = values_tmp(2);
            J_values_final(end + 1) = values_tmp(3);
            dS_values_final(end + 1) = values_tmp(4);
            dA_values_final(end + 1) = values_tmp(5);
        end % while loop
    end % if variables don't exist
    
    % Post-process IDS and IDA to discover number of rows & columns in
    % original parameter sweep data
    initial_data_parameter_s_values = unique(IDS);
    n_initial_data_parameter_s_values = length(initial_data_parameter_s_values);
    initial_data_parameter_a_values = unique(IDA);
    n_initial_data_parameter_a_values = length(initial_data_parameter_a_values);
    
    IDS = reshape(IDS, n_initial_data_parameter_s_values, n_initial_data_parameter_a_values);
    IDA = reshape(IDA, size(IDS));
    J_values_final = reshape(J_values_final, size(IDS));
    dS_values_final = reshape(dS_values_final, size(IDS));
    dA_values_final = reshape(dA_values_final, size(IDS));
    
    %% Create 2D plots
    figure('Name', '2D Initial Declination Plots')
    % J(v_{final})
    subplot(3,1,1);
    imagesc( ...
        'XData', initial_data_parameter_a_values, ...
        'YData', initial_data_parameter_s_values, ...
        'CData', J_values_final);
    xlabel('\lambda_a');
    ylabel('\lambda_s');
    title('Initial Declination vs. J(v_{final})');
    colorbar
    
    % ||s_{final} - s_{true}||
    subplot(3,1,2)
    imagesc( ...
        'XData', initial_data_parameter_a_values, ...
        'YData', initial_data_parameter_s_values, ...
        'CData', dS_values_final);
    xlabel('\lambda_a');
    ylabel('\lambda_s');
    title('Initial Declination vs. ||s_{final} - s_{true}||');
    colorbar
    
    % ||a_{final} - a_{true}||
    subplot(3,1,3)
    imagesc( ...
        'XData', initial_data_parameter_a_values, ...
        'YData', initial_data_parameter_s_values, ...
        'CData', dA_values_final);
    xlabel('\lambda_a');
    ylabel('\lambda_s');
    title('Initial Declination vs. ||a_{final} - a_{true}');
    colorbar
    
    % Add overall plot title
    if exist('sgtitle')
      sgtitle(sprintf('Nx=%d,Nt=%d', len_xmesh, len_tmesh));
    else
      title(sprintf('Nx=%d,Nt=%d', len_xmesh, len_tmesh));
    end
    saveas(gcf, plot_filename_2d, 'pdf');
    
    % 1D (slice) plots
    figure('Name', '1D Plot Examples');
    subplot(3,1,1);
    % J(v_{final}) for selected lambda_a (so a fixed column)
    initial_data_parameter_a = IDA(1,1);
    plot(IDS(:, 1), J_values_final(:, 1));
    xlabel('\lambda_s')
    ylabel('J(v_{final})')
    title(sprintf('Dependence of J(v_{final}) on \\lambda_s for \\lambda_a == %2.5f', initial_data_parameter_a));
    
    % ||s_{final} - s_{true}|| for selected lambda_a (so a fixed column)
    subplot(3,1,2);
    initial_data_parameter_a = IDA(1,1);
    plot(IDS(:, 1), dS_values_final(:, 1));
    xlabel('\lambda_s')
    ylabel('||s_{final} - s_{true}||')
    title(sprintf('Dependence of ||s_{final} - s_{true}|| on \\lambda_s for \\lambda_a == %2.5f', initial_data_parameter_a));
    
    % ||a_{final} - a_{true}|| for selected lambda_a (so a fixed column)
    subplot(3,1,3);
    initial_data_parameter_a = IDA(1,1);
    plot(IDS(:, 1), dA_values_final(:, 1));
    xlabel('\lambda_s')
    ylabel('||a_{final} - a_{true}||')
    title(sprintf('Dependence of ||a_{final} - a_{true}|| on \\lambda_s for \\lambda_a == %2.5f', initial_data_parameter_a));
    
    % Add overall plot title
    if exist('sgtitle')
      sgtitle(sprintf('Nx=%d,Nt=%d', len_xmesh, len_tmesh));
    else
      title(sprintf('Nx=%d,Nt=%d', len_xmesh, len_tmesh));
    end
    saveas(gcf, plot_filename_1d, 'pdf');
end % if generate_plots
end
