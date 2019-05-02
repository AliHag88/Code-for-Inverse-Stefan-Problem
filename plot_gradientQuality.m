function [] = plot_gradientQuality(generate_data, generate_plots, Nx, Nt, n_fd_epsilon)
  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.
  optimization_parameter_defaults;

  % Set default values of the input arguments
  if ~exist('generate_data', 'var')
    generate_data = true;
  end
  if ~exist('generate_plots', 'var')
    generate_plots = false;
  end
                                % Number of iterations to test
  if ~exist('n_fd_epsilon', 'var')
    n_fd_epsilon = 20;
  end

  % Compute filename corresponding to this example
  output_filename = sprintf('kappa_Nx%d_Nt%d.csv', Nx, Nt);
  plot_filename = sprintf('kappa_Nx%d_Nt%d.pdf', Nx, Nt);

  % Headers and file format for output
  output_headers = ['fd_epsilon,kappa_s,kappa_s_precond,kappa_a,kappa_a_precond,' ...
                      'logkappa_s,logkappa_s_precond,logkappa_a,logkappa_a_precond\n' ];

  % Output format and size of each row
  output_format = '%.16f,';
  output_size = [1, 9];

  % Generate data
  if generate_data

    % Initialize storage
    fd_epsilon_values = linspace(eps(0), 1e-1, n_fd_epsilon)';
    kappa_s = zeros(size(fd_epsilon_values))*NaN;
    kappa_s_precond = zeros(size(fd_epsilon_values))*NaN;
    kappa_a = zeros(size(fd_epsilon_values))*NaN;
    kappa_a_precond = zeros(size(fd_epsilon_values))*NaN;

    % Open file to export data in CSV format
    fileID = fopen(output_filename, 'w');

    % Add header
    fprintf(fileID, output_headers);

    % Calculate kappa values
    for i = 1:length(fd_epsilon_values)
      fd_epsilon_curr = fd_epsilon_values(i);
      try
        [kappa_s(i), kappa_s_precond(i), kappa_a(i), kappa_a_precond(i)] = ...
        gradientQuality(Nx, Nt, fd_epsilon_curr, ...
                        initial_data_parameter_s, initial_data_parameter_a, ...
                        regularization_s, regularization_a);
      catch ME # Fail with grace (by ignoring the issue.)
        disp(ME)
      end

      % Output to file
      fprintf(fileID, [repmat(output_format, output_size) '\n'], ...
              fd_epsilon_values(i), ...
              kappa_s(i), kappa_s_precond(i), ...
              kappa_a(i), kappa_a_precond(i), ...
              log10(abs(kappa_s(i) - 1)), log10(abs(kappa_s_precond(i) - 1)), ...
              log10(abs(kappa_a(i) - 1)), log10(abs(kappa_a_precond(i) - 1)) ...
             );
    end
    % Close data file
    fclose(fileID);

  end % if generate_data


  % Generate plots
  if generate_plots
  % If data already exists in this workspace, use it.
  % Otherwise, read it from the file.
    if ~exist('kappa_s', 'var') && isfile(output_filename)
      fileID = fopen(output_filename, 'r');
      fgetl(fileID); % Skip header

      % Initialize storage
      fd_epsilon_values = [];
      kappa_s = [];
      kappa_s_precond = [];
      kappa_a = [];
      kappa_a_precond = [];

      % Read all lines from file
      while true
        tline = fgetl(fileID);
        if tline == -1 % Documented as EOL marker from fgetl
          break
        end

        % Parse data from current line
        values_tmp = sscanf(tline, '%f,', output_size);

        % Extract values from temp array
        fd_epsilon_values(end + 1) = values_tmp(1);
        kappa_s(end + 1) = values_tmp(2);
        kappa_s_precond(end + 1) = values_tmp(3);
        kappa_a(end + 1) = values_tmp(4);
        kappa_a_precond(end + 1) = values_tmp(5);
      end % while loop
    end % if variables don't exist

    %% Plots for s(t)
    % kappa_s
    subplot(2,2,1);
    plot(fd_epsilon_values, kappa_s);
    xlabel('\epsilon');
    ylabel('\kappa_s');
    title(sprintf('\\kappa_s, Precond:False', Nx, Nt));

    % kappa_s_precond
    subplot(2,2,2);
    plot(fd_epsilon_values, kappa_s_precond);
    xlabel('\epsilon');
    ylabel('\kappa_s');
    title(sprintf('\\kappa_s, Precond:True', Nx, Nt));

    %% Plots for a(t)
    % kappa_a
    subplot(2,2,3);
    plot(fd_epsilon_values, kappa_a);
    xlabel('\epsilon');
    ylabel('\kappa_a');
    title(sprintf('\\kappa_a, Precond:False', Nx, Nt));

    % kappa_a_precond
    subplot(2,2,4);
    plot(fd_epsilon_values, kappa_a_precond);
    xlabel('\epsilon');
    ylabel('\kappa_a');
    title(sprintf('\\kappa_a, Precond:True', Nx, Nt));

    if exist('sgtitle')
      sgtitle(sprintf('Nx=%d,Nt=%d', Nx, Nt));
    else
      title(sprintf('Nx=%d,Nt=%d', Nx, Nt));
    end
    saveas(gcf, plot_filename, 'pdf');
  end % if generate_plots

end % Function
