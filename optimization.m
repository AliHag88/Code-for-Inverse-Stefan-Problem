function [J, svals, avals] = optimization(len_xmesh, len_tmesh, tolerance, num_iterations, do_visualization)
  % optimization: Run ISP example
  % Input Arguments:
  %    - len_xmesh: Number of space grid points. Default: 20
  %    - len_tmesh: Number of time grid points. Default: 20
  %    - tolerance: Required error for stopping criteria. Default: 1e-5
  %    - num_iterations: Number of gradient descent steps to take. Default: 200
  %    - do_visualization: Set to true to emit visualizations during the
  %      optimization process.

  % Set default arguments
  if ~exist('len_xmesh', 'var')
    len_xmesh = 50;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 10;
  end
  if ~exist('tolerance', 'var')
    tolerance = 1e-5;
  end
  if ~exist('num_iterations', 'var')
    num_iterations = 500;
  end
  if ~exist('do_visualization', 'var')
    do_visualization = true;
  end

  % Set final moment (eliminate one "magic constant" used in subsequent locations.)
  t_final = 1;

  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0,1,len_xmesh);  % Space discretization (row)
  tmesh = linspace(0,t_final,len_tmesh)'; % Time discretization (column)

  % Initial setup for solver (all tunable parameters should be set here)
  [step_size, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh);

  % Initialize svals and avals
  svals = s_ini(tmesh);
  s_der = est_deriv(svals, tmesh);
  avals = a_ini(tmesh);

  % Vector of cost functional values
  J = zeros(num_iterations,1);

  % Counter for number of iterations
  k = 0;

  % Main Optimization Loop
  while k < num_iterations
    k = k + 1;

    % Calculate solution of forward problem
    [au_xx_S, u_x_S, u_S, u_T, u] = ...
      Forward(xmesh, tmesh, svals, avals, g, u_true_0);

    % Calculate solution of adjoint problem
    [psi_t_S, psi_x_S, psi_S, psi_T, psi] = ...
      Adjoint(xmesh, tmesh, svals, avals, u_T, w_meas, s_der, u_S, mu_meas);

    % Take step in antigradient direction
    s_update = grad_s(tmesh, s_der, svals, w_meas, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_S, u_S, au_xx_S, s_star, psi_T);
    s_update = s_update / norm(s_update);

    svals = svals - step_size * s_update;

    avals = avals - 0 * grad_a(u, psi, tmesh); % Note: avals not updated.

    % Precondition gradient
    % grad=precond(tmesh, grad, L); % Preconditioning for s(t) gradient

    % Update s_der
    s_der = est_deriv(svals, tmesh);

    % Calculate functional
    J(k) = abs(svals(end) - s_star)^2 + ...
           trapz(xmesh, (u_T - w_meas(xmesh)).^2) + ...
           trapz(tmesh, (u_S - mu_meas(tmesh)).^2);

    % Do visualization if selected
    if do_visualization
        pause_time = 0; % Second
        visualization(xmesh, tmesh, svals, avals, u, k, J, pause_time);
        drawnow();
    end

    % Check stopping criteria
    if J(k) < tolerance
        J = J(1:k);
        break
    end

  end % End main loop

end % optimization function
