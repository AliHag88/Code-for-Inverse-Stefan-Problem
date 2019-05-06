function [] = plot_a_and_s_with_precond(len_xmesh, len_tmesh, n_sobolev_preconditioning_s_values,...
                                   n_sobolev_preconditioning_a_values,tolerance, num_iterations,...
                                   num_sub_iterations, use_synthetic_data,initial_data_parameter_s,...
                                   initial_data_parameter_a,n_initial_data_parameter_s_values,...
                                   n_initial_data_parameter_a_values)
  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.)
  optimization_parameter_defaults;

  % Set more default values for parameters
if ~exist('initial_data_parameter_s', 'var')
  n_initial_data_parameter_s_values = 100;
end
if ~exist('initial_data_parameter_a', 'var')
  n_initial_data_parameter_a_values = 100;
end


  % Endpoints of Convergence Interval for
  % initial_data_parameter_s/initial_data_parameter_a parameters
 
  a1=-1;  % left endpoint or l_a search
  a2=1;  % right endpoint or l_a search
  
  s1=-1;  % left endpoint for l_s search
  s2=1;  % right endpoint for l_s search
   


    % Choose range of sobolev_preconditioning_s values
    initial_data_parameter_a_values = linspace(a1, a2, n_initial_data_parameter_a_values)';

    initial_data_parameter_s_values = linspace(s1, s2, n_initial_data_parameter_s_values)';

    
    % Initialize storage for output data
    J_values_final = zeros(n_initial_data_parameter_a_values+1,n_sobolev_preconditioning_s_values+1);
    dAS_values_final = zeros(n_initial_data_parameter_a_values+1,n_sobolev_preconditioning_s_values+1);
 

    % Discretization for both for forward and adjoint problem
    tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

    % Calculate corresponding s_true and a_true
    [~, ~, ~, ~, ~, s_true, ~, a_true] = true_solution(tmesh);

    % Since the mesh is fixed, we can compute the analytic solution once.
    s_true_values = s_true(tmesh)';
    a_true_values = a_true(tmesh)';

 
    
    for i = 2:n_initial_data_parameter_a_values+1
      initial_data_parameter_a = initial_data_parameter_a_values(i-1);
       for j = 2:n_initial_data_parameter_s_values+1
         initial_data_parameter_s = initial_data_parameter_s_values(i-1);
     
        [jvals, svals, avals] = ...
          optimization( ...
                        len_xmesh, len_tmesh, ...
                        tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
                        initial_data_parameter_s, initial_data_parameter_a, ...
                        sobolev_preconditioning_s, sobolev_preconditioning_a, ...
                        reconstruct_s, reconstruct_a, ...
                        false);

 J_values_final(i,j) = jvals(end);
 dAS_values_final(i,j) = sqrt(trapz(tmesh, (svals(end, :) - s_true_values).^2))/sqrt(trapz(tmesh, (s_true_values).^2))+...
 sqrt(trapz(tmesh, (avals(end, :) - a_true_values).^2))/sqrt(trapz(tmesh, (a_true_values).^2));   
     end
    end  % Loop over sobolev_preconditioning_s_values
    
    % adding l_s and l_a vectors to the matrix 
    J_values_final(1,2:end)=initial_data_parameter_s_values;
    J_values_final(2:end,1)=initial_data_parameter_a_values;
    
    dAS_values_final(1,2:end)=initial_data_parameter_s_values;
    dAS_values_final(2:end,1)=initial_data_parameter_a_values;    
    
           
    % Transferring matrix values to the txt files   
    dlmwrite('J_values.txt',J_values_final,'delimiter', '\t','newline', 'pc') 
    dlmwrite('dSA_values.txt',dAS_values_final,'delimiter', '\t','newline', 'pc') 
    
    % Printing elements of l_s and l_a
    fprintf('row elements correspond to lambda_s = %s\n', sprintf('%d ', initial_data_parameter_s_values))
    fprintf('column elements correspond to lambda_a = %s\n', sprintf('%d ', initial_data_parameter_a_values))
    
    
% Ploting 2D convergence analysis for a&s

imagesc('XData',initial_data_parameter_s_values,'YData',initial_data_parameter_a_values,'CData',J_values_final(2:end,2:end))
colorbar
axis image
caxis([0 5])

   
imagesc('XData',initial_data_parameter_s_values,'YData',initial_data_parameter_a_values,'CData',dAS_values_final(2:end,2:end))
colorbar
axis image
caxis([0 5])



end