function [estimate_s_precond, estimate_s, estimate_a_precond, estimate_a] = ...
         gradientQuality(len_xmesh, len_tmesh, ...
                         fd_epsilon, ...
                         initial_data_parameter_s, initial_data_parameter_a, ...
                         sobolev_preconditioning_s, sobolev_preconditioning_a)
  % Calculate an estimate of the gradient quality by comparing our analytic
  % gradient to values estimated from a finite difference about s_ini and a_ini
  % as created by initial_setup, in an arbitrarily chosen direction.

  % Call same setup code as in optimization script
  % (any duplicated variables will be overwritten later.
  optimization_parameter_defaults;
  
  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0, 1, len_xmesh);  % Space discretization (row)
  tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)

  % Initial setup for solver (all tunable parameters should be set here)
  [~, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh, use_synthetic_data, initial_data_parameter_s, initial_data_parameter_a);

  % Initialize svals and avals
  s_old = s_ini(tmesh);
  a_old = a_ini(tmesh);

  % Create search direction for s and a
  dv = ones(size(tmesh));
  dv(1) = 0;

  % Calculate solution of forward problem and functional value at initial approach
  [au_xx_S, u_x_S, u_S, u_T, u, J_value] = ...
    Functional(xmesh, tmesh, s_old, a_old, g, u_true_0, s_star, w_meas, mu_meas);

  % Calculate solution of adjoint problem
  [psi_t_S, psi_x_S, psi_S, psi] = ...
    Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);

  % Calculate update direction vector s_update
  [s_update_1, s_update_1_T] = grad_s(tmesh, s_old, w_meas, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_S, u_S, au_xx_S, s_star);

  a_update_1 = grad_a(u,psi,tmesh,xmesh,s_old);

  % Preconditioning for s(t) gradient
  s_update_2 = precond( ...
                        s_precond_mode, sobolev_preconditioning_s, tmesh, ...
                        s_update_1, s_update_1_T ...
                      );

  % Preconditioning for a(t) gradient
  a_update_2 = precond( ...
                        a_precond_mode, sobolev_preconditioning_a, tmesh, ...
                        a_update_1 ...
                      );

  % Estimate finite difference using s_old
  [~, ~, ~, ~, ~, J_value_s_plus] = ...
    Functional(xmesh, tmesh, s_old + fd_epsilon * dv, a_old, g, u_true_0, s_star, w_meas, mu_meas);
  fd_est_s = (J_value_s_plus - J_value) / fd_epsilon;

  % Estimate finite difference using a_old
  [~, ~, ~, ~, ~, J_value_a_plus] = ...
    Functional(xmesh, tmesh, s_old, a_old + fd_epsilon * dv, g, u_true_0, s_star, w_meas, mu_meas);
  fd_est_a = (J_value_a_plus - J_value) / fd_epsilon;

  % Calculate estimates for ratio of directional gradients
  estimate_s = fd_est_s / dot(s_update_1, dv);
  estimate_a = fd_est_a / dot(a_update_1, dv);

  % Calculate estimates for ratio of directional gradients
  estimate_s_precond = fd_est_s / dot(s_update_2, dv);
  estimate_a_precond = fd_est_a / dot(a_update_2, dv);

  fprintf('fd_epsilon=%5.10f, kappa_s=%5.10f, kappa_s_precond=%5.10f, kappa_a=%5.10f, kappa_a_precond=%5.10f\n', ...
          fd_epsilon, estimate_s, estimate_s_precond, estimate_a, estimate_a_precond);
end
