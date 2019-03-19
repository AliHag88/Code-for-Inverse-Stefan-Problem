function [J_values, s_values, a_values] = optimization(len_xmesh, len_tmesh, tolerance, num_iterations, num_sub_iterations, do_visualization)
  % optimization: Run ISP example
  % Input Arguments:
  %    - len_xmesh: Number of space grid points. Default: 20
  %    - len_tmesh: Number of time grid points. Default: 20
  %    - tolerance: Required error for stopping criteria. Default: 1e-5
  %    - num_iterations: Number of gradient descent steps to take. Default: 200
  %    - num_sub_iterations: Number of trial steps to take to find decreasing step. Default: 20
  %    - do_visualization: Set to true to emit visualizations during the
  %      optimization process.
  % Output Arguments:
  %    - J_values: Vector of functional values at each gradient iteration.
  %    - s_values: Vector of controls x=s_k(t) at each gradient iteration.
  %    - a_values: Vector of controls a_k(t) at each gradient iteration.
  % Note that for the last two output arguments, columns correspond to the iteration number
  % and each row is a separate control vector.

  % Set default arguments
  if ~exist('len_xmesh', 'var')
    len_xmesh = 20;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 20;
  end
  if ~exist('tolerance', 'var')
    tolerance = 1e-5;
  end
  if ~exist('num_iterations', 'var')
    num_iterations = 200;
  end
  if ~exist('num_sub_iterations', 'var')
    num_sub_iterations = 20;
  end
  if ~exist('do_visualization', 'var')
    do_visualization = false;
  end
  
  % If the norm of the update vector is below the threshold below, we will not normalize it.
  norm_update_threshold = 1e-10;
  
  % Thresholds for minimum values of s(t) and a(t)
  svals_minimum_threshold = 1e-4;
  avals_minimum_threshold = 1e-4;
  
  % Counter for number of iterations
  k = 1;
  
  % Set final moment
  t_final = 1;
  
  % Vector of cost functional values
  J_values = zeros(num_iterations + 1, 1);
  
  % Vector of boundary curve iterates
  s_values = zeros(num_iterations + 1, len_tmesh);
  
  % Vector of diffusion coefficient iterates;
  a_values = zeros(num_iterations + 1, len_tmesh);

  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0, 1, len_xmesh);  % Space discretization (row)
  tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

  % Initial setup for solver (all tunable parameters should be set here)
  [max_step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh);

  % Initialize svals and avals
  s_old = s_ini(tmesh);
  a_old = a_ini(tmesh);

  % Store initial approach to s_values and a_values
  s_values(k, :) = s_old;
  a_values(k, :) = a_old;
  
  % Calculate solution of forward problem and functional value at initial approach
  [au_xx_S, u_x_S, u_S, u_T, u, J_values(k)] = ...
    Functional(xmesh, tmesh, s_old, a_old, g, u_true_0, s_star, w_meas, mu_meas);

  % Calculate solution of adjoint problem
  [psi_t_S, psi_x_S, psi_S, psi_T, ~] = ...
    Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);
  
  
  % Main Optimization Loop
  while k <= num_iterations
    k = k + 1;
    if any(isnan(psi_t_S)) || any(isnan(psi_x_S)) || any(isnan(psi_S)) || any(isnan(psi_T)) ...
          || any(isnan(au_xx_S)) || any(isnan(u_x_S)) || any(isnan(u_S)) || any(isnan(u_T))
      disp(['Invalid state or adjoint at step ' num2str(k) '.']);
      break
    end
      
    % Calculate update direction vector s_update
    s_update = grad_s(tmesh, s_old, w_meas, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_S, u_S, au_xx_S, s_star, psi_T);
    if norm(s_update) > norm_update_threshold % Or is generally too small to work with
        s_update = s_update / norm(s_update);
    end
    
    % Calculate update direction vector a_update
    %a_update = grad_a(u, psi, tmesh);
    %if norm(a_update) > 1e-10
    %    a_update = a_update / norm(a_update);
    %end

    % Preconditioning for s(t) gradient
    % s_update = precond(tmesh, s_update, L);
    
    % Preconditioning for a(t) gradient
    % a_update = precond(tmesh, a_update, L);    
    
    curr_step_size = max_step_size;
    sub_iter = 1;
    while true
      % Take trial step along direction vector s_update and a_update
      s_new = s_old - curr_step_size * s_update;
  
      a_new = a_old;% - 0 * a_update; % Note: avals not updated.
      
      % If svals or avals becomes to small, reduce the step size and try again
      if any(s_new < svals_minimum_threshold) || any(a_new < avals_minimum_threshold)  
        if sub_iter >= num_sub_iterations
            disp(['Singularity developed in s(t) or a(t) at iteration ' num2str(sub_iter) '. Stopping Iterative Method.']);
            disp('s_final: ');
            disp(s_new)
            break
        end
        disp(['Reducing trial step size to ' num2str(curr_step_size/2)]);
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
        disp(['Found a decreasing step after ' num2str(sub_iter) ' iterations.']);
        disp(['Functional value J_{' num2str(k) '} == ' num2str(J_curr) '.']);
        J_values(k) = J_curr;
        s_old = s_new; s_values(k, :) = s_old;
        a_old = a_new; a_values(k, :) = a_old;
        break
      end
      
      % If we can't find a step size that decreases the functional value,
      if sub_iter >= num_sub_iterations
        disp(['Failed to find appropriate step size in ' num2str(sub_iter) ' iterations.']);
        disp('s_final: ');
        disp(s_new)
        return
      end
      
      % Reduce step size and try again
      curr_step_size = curr_step_size / 2;
      sub_iter = sub_iter + 1;
    end % While loop for sub step
    
    % Do visualization if selected
    if do_visualization
        pause_time = 0; % Units for pause_time are seconds
        visualization(xmesh, tmesh, s_old, a_old, u, k, J_values, pause_time);
        drawnow();
    end

    % Check stopping criteria
    if (J_values(k) - J_values(k-1)) < tolerance * J_values(k)
        disp(['Iterations stationary (in relative error) at k=' num2str(k) ' with tolerance ' num2str(tolerance) ]);
        J_values = J_values(1:k);
        break
    end

    % Update adjoint
    [psi_t_S, psi_x_S, psi_S, psi_T, ~] = ...
        Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);
  end % End main loop

end % optimization function
