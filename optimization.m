function [J, svals, avals] = optimization(len_xmesh, len_tmesh, tolerance, num_iterations, do_visualization)
  % optimization: Run ISP example
  % Input Arguments:
  %    - len_xmesh: Number of space grid points. Default: 100
  %    - len_tmesh: Number of time grid points. Default: 100
  %    - tolerance: Required error for stopping criteria. Default: 1e-5
  %    - num_iterations: Number of gradient descent steps to take. Default: 200
  %    - do_visualization: Set to true to emit visualizations during the
  %      optimization process.

  % Set default arguments
  if ~exist('len_xmesh', 'var')
    len_xmesh = 100;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 100;
  end
  if ~exist('tolerance', 'var')
    tolerance = 1e-5
  end
  if ~exist('num_iterations', 'var')
    num_iterations = 200
  end
  if ~exist('do_visualization', 'var')
    do_visualization = false;
  end

  % Set final moment (eliminate one "magic constant" used in subsequent locations.)
  t_final = 1;

  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0,1,len_xmesh);  % Space discretization
  tmesh = linspace(0,t_final,len_tmesh); % Time discretization

  % Initial setup for solver (contains all parameters)
  % Uses global information to set mu_meas, w_meas, alpha
  % as well as initial approach s_ini and a_ini
  %L=1; % Preconditioning 
  initial_setup;

  % Initialize svals and avals
  svals = s_ini(tmesh);
  s_der = deriv_est(svals, tmesh);
  avals = a_ini(tmesh);

  % Vector of cost functional values
  J = zeros(num_iterations,1);

  % Counter for number of iterations
  k = 0;

  % Main Optimization Loop
  while k < num_iterations
    k = k + 1

    % Calculate solution of forward problem
    [au_xx_S, u_x_S, u_S, u_T, u] = ...
      Forward(xmesh, tmesh, svals, avals, u_true_0);

    % Calculate solution of adjoint problem
    [psi_T, psi_t, psi_t_S, psi_S, psi_x_S, psi_x, psi] = ...
      Adjoint(xmesh, tmesh, svals, avals, u_T, w_meas, s_der, u_S, mu_meas);

    % Take step in antigradient direction
    svals = svals - alpha * grad_s(u, psi, psi_x, tmesh, s_der, svals, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_t, psi_S, u_S, au_xx_S, s_star, psi_T);
    avals = avals - 0 * grad_a(u,psi,tmesh); % Note: avals not updated.


    % Precondition gradient
    % grad=precond()'; % Preconditioning for s(t) gradient

    % Update s_der
    s_der = deriv_est(svals);

    % Calculate functional
    J(k) = abs(svals(end) - s_star)^2 + ...
           trapz(xmesh, (u_T - w_meas(xmesh)).^2) + ...
           trapz(tmesh, (u_S - mu_meas(tmesh)).^2);

    % Do visualization if selected
    if do_visualization
      visualization(xmesh, tmesh, svals, avals, u, u_true, s_true, k, J)
    end

    # Check stopping criteria
    if J(k) < tolerance
      break
    end

  end

end # optimization function
