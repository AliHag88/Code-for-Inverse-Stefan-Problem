function [J_values, s_values, a_values] = optimization(...
    len_xmesh, len_tmesh, ...
    tolerance, num_iterations, num_sub_iterations, use_synthetic_data, ...
    initial_data_parameter_s, initial_data_parameter_a, ...
    sobolev_preconditioning_s, sobolev_preconditioning_a, ...
    reconstruct_s, reconstruct_a, do_visualization)
  % optimization: Run ISP example
  % Input Arguments:
  %    - len_xmesh: Number of space grid points. Default: 20
  %    - len_tmesh: Number of time grid points. Default: 20
  %    - tolerance: Required error for stopping criteria. Default: 1e-5
  %    - num_iterations: Number of gradient descent steps to take. Default: 200
  %    - num_sub_iterations: Number of trial steps to take to find decreasing step. Default: 20
  %    - use_synthetic_data: Set true to use measurements from PDE solver,
  %      or set to false to use measurements from given problem. Default: true
  %    - initial_data_parameter_{s,a}: Parameter passed into initial_setup
  %      to control how far initial approach is from analytic solution.
  %      Default: 0 (s_initial == s_true or a_initial == a_true)
  %    - regularization_{s,a}: Weight given to regularization term during
  %      gradient descent process.
  %    - reconstruct_{s,a}: Boolean values to select whether to reconstruct s
  %      and/or a. If false, the initial approach will be used at each iteration.
  %    - do_visualization: Set to true to emit visualizations during the optimization process.
  %      Default: false
  % Output Arguments:
  %    - J_values: Vector of functional values at each gradient iteration.
  %    - s_values: Vector of controls x=s_k(t) at each gradient iteration.
  %    - a_values: Vector of controls a_k(t) at each gradient iteration.
  % Note that for the last two output arguments, columns correspond to the iteration number
  % and each row is a separate control vector.

  % Set default arguments
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
    num_sub_iterations = 6;
  end
  if ~exist('do_visualization', 'var')
    do_visualization = true;
  end
  
  % See initial_setup.m for the following parameters
  if ~exist('use_synthetic_data', 'var')
    use_synthetic_data = true;
  end
  if ~exist('initial_data_parameter_s', 'var')
    initial_data_parameter_s = 1.9;
  end
  if ~exist('initial_data_parameter_a', 'var')
    initial_data_parameter_a = 0;
  end
  % Preconditioning parameters
  if ~exist('sobolev_preconditioning_s', 'var')
    sobolev_preconditioning_s = 0.5;
  end
  if ~exist('sobolev_preconditioning_a', 'var')
    sobolev_preconditioning_a = 0.22;
  end
  % Choosing coefficients for reconstruction
  if ~exist('reconstruct_a', 'var')
    reconstruct_a = 1;
  end
   if ~exist('reconstruct_s', 'var')
    reconstruct_s = 1;
  end
  
  
  % If the norm of the update vector is below the threshold below, we will not normalize it.
  norm_update_threshold = 1e-10;
  
  % Thresholds for minimum values of s(t) and a(t)
  svals_minimum_threshold = 1e-4;
  avals_minimum_threshold = 1e-4;

  % Preconditioning "mode".
  % Select == 1 to impose homogeneous Neumann condition (the "original" flavor)
  % Select == 2 to give a Neumann condition on the right-hand side and a
  %             homogeneous Dirichlet condition on the left-hand side
  s_precond_mode = 2;
  a_precond_mode = 2;
  
  % Counter for number of iterations
  k = 1;
  
  % Set final moment
  t_final = 1;
  
  % Vector of cost functional values. NaNs are used as a
  % placeholder.
  J_values = zeros(num_iterations + 1, 1) * NaN;
  
  % Vector of boundary curve iterates
  s_values = zeros(num_iterations + 1, len_tmesh) * NaN;
  
  % Vector of diffusion coefficient iterates;
  a_values = zeros(num_iterations + 1, len_tmesh) * NaN;

  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0, 1, len_xmesh);  % Space discretization (row)
  tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

  % Take "true" (analytic) values of s and a
  [~, ~, ~, ~, ~, s_true, ~, a_true] = true_solution(tmesh);

  % Calculate values of analytic solution on time grid
  s_true_values = s_true(tmesh);
  a_true_values = a_true(tmesh);
  
  % Initial setup for solver (all tunable parameters should be set here)
  [max_step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh, use_synthetic_data, initial_data_parameter_s, initial_data_parameter_a);

  % Initialize svals and avals
  s_old = s_ini(tmesh);
  a_old = a_ini(tmesh);

  % Store initial approach to s_values and a_values
  s_values(k, :) = s_old;
  a_values(k, :) = a_old;
  
  % Calculate solution of forward problem and functional value at initial approach
  [au_xx_S, u_x_S, u_S, u_T, u, J_values(k)] = ...
    Functional(xmesh, tmesh, s_old, a_old, g, u_true_0, s_star, w_meas, mu_meas);

  fprintf('Initial functional value: %2.5f.\n', J_values(k));

  % Calculate solution of adjoint problem
  [psi_t_S, psi_x_S, psi_S, psi] = ...
    Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);
  
  
  % Main Optimization Loop
  while k <= num_iterations
    k = k + 1;
    if any(isnan(psi_t_S)) || any(isnan(psi_x_S)) || any(isnan(psi_S)) ...
          || any(isnan(au_xx_S)) || any(isnan(u_x_S)) || any(isnan(u_S))
      fprintf('Invalid state or adjoint at step %d.\n', k);
      break
    end
      
    % Calculate update direction vector s_update
    s_update, s_update_T = grad_s(tmesh, s_old, w_meas, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_S, u_S, au_xx_S, s_star);
    % Only normalize if the update vector has numerically nonzero norm.
    if norm(s_update) > norm_update_threshold
        s_update = s_update / norm(s_update);
    end
    
    % Calculate update direction vector a_update
    a_update = grad_a(u,psi,tmesh,xmesh,s_old);
    if norm(a_update) > 1e-10
        a_update = a_update / norm(a_update);
    end
   
    % Preconditioning for s(t) gradient
    s_update = precond(s_precond_mode, sobolev_preconditioning_s, tmesh, s_update, s_update_T);
    
    % Preconditioning for a(t) gradient
    a_update = precond(a_precond_mode, sobolev_preconditioning_a, tmesh, a_update);
    
    curr_step_size = max_step_size;
    sub_iter = 1;
    while true
      % Take trial step along direction vector s_update and a_update
      s_new = s_old - reconstruct_s * curr_step_size * s_update;
  
      a_new = a_old - reconstruct_a * 0.01 * a_update; % Note: avals not updated.
                                     % Note: the above comment is obviously
                                     % incorrect.
      
      % If svals or avals becomes to small, reduce the step size and try again
      if any(s_new < svals_minimum_threshold) || any(a_new < avals_minimum_threshold)  
        if sub_iter >= num_sub_iterations
            fprintf('Singularity developed in s(t) or a(t) after %d sub-iteration(s).\n', sub_iter);
            fprintf('Stopping gradient descent at k=%d.\n', k);
            disp('s_final: ');
            disp(s_new)
            break
        end
        curr_step_size = curr_step_size / 2;
        sub_iter = sub_iter + 1;
        continue
      end
      
      % Calculate functional value and state vector at new svals and avals vectors
      [au_xx_S, u_x_S, u_S, u_T, u, J_curr] = ...
        Functional(xmesh, tmesh, s_new, a_new, ...
                   g, u_true_0,  s_star, w_meas, mu_meas...
                   );

      % If we've found a step that decreases the functional value, break out of
      % the loop after saving s_new and a_new over s_old and a_old.
      if J_curr < J_values(k-1)
        fprintf('Found a decreasing step after %d sub-iterations.\n', sub_iter);
        fprintf('Functional value J_{%d} == %2.5f.\n', k, J_curr);
        fprintf('||s_k-s_true||/||s_true||=%2.5f.\n', norm(s_new-s_true_values)/norm(s_true_values));
        fprintf('||a_k-a_true||/||a_true||=%2.5f.\n', norm(a_new-a_true_values)/norm(a_true_values));

        J_values(k) = J_curr;
        s_old = s_new; s_values(k, :) = s_old; s_values = s_values(1:k, :);
        a_old = a_new; a_values(k, :) = a_old; a_values = a_values(1:k, :);
        break
      end
      
      % If we can't find a step size that decreases the functional value,
      % Save the final iterate and bail out of the gradient descent process.
      if sub_iter >= num_sub_iterations
        fprintf('Failed to find appropriate step size in %d sub-iterations.\n', sub_iter);
        fprintf('Functional value J_{%d} == %2.5f.\n', k, J_curr);
        fprintf('||s_k-s_true||/||s_true||=%2.5f.\n', norm(s_new-s_true_values)/norm(s_true_values));
        fprintf('||a_k-a_true||/||a_true||=%2.5f.\n', norm(a_new-a_true_values)/norm(a_true_values));
        disp('s_final: ');
        disp(s_new);

        s_values(k, :) = s_new; s_values = s_values(1:k, :);
        a_values(k, :) = a_new; a_values = a_values(1:k, :);
        J_values = J_values(1:k);
        return
      end
      
      % Reduce step size and try again. This step size selection method is
      % precisely the Armijo rule with bisection at each trial step.
      curr_step_size = curr_step_size / 2;
      sub_iter = sub_iter + 1;
    end % While loop for sub step
    
    % Do visualization if selected
    if do_visualization
        pause_time = 0; % Units for pause_time are seconds
        visualization(xmesh, tmesh, s_old, a_old, u, k, J_values, pause_time, initial_data_parameter_s, initial_data_parameter_a);
        drawnow();
    end

    % Check stopping criteria
    if abs(J_values(k) - J_values(k-1)) < tolerance * J_values(k)% && abs(a_new - a_old) < tolerance * a_new
        fprintf('Iterations stationary (in relative error) at k=%d with tolerance %2.5f.\n', k, tolerance);
        J_values = J_values(1:k);
        break
    end

    % Update adjoint
    [psi_t_S, psi_x_S, psi_S, ~] = ...
        Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);
  end % End main loop

end % optimization function
