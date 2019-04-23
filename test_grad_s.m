function [error_out_grad_s, error_out_precond_s] = test_grad_s(len_xmesh, len_tmesh, generate_validation)
%TEST_GRAD_S Runs integration test for grad_s and precond_s codes
% for now, actually checks for norm of declination between calculated and
% stored gradient and preconditioned gradient
  
  if ~exist('len_xmesh', 'var')
    len_xmesh = 10;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 10;
  end
  if ~exist('generate_validation', 'var')
    generate_validation = false;
  end
  
  t_final = 1;
  use_synthetic_data = true;
  initial_data_parameter_s = 0.6;
  initial_data_parameter_a = 0.6;
  norm_update_threshold = 1e-10;
  L_s = 0.6;
  
  % Discretization for both for forward and adjoint problem
  xmesh = linspace(0, 1, len_xmesh);  % Space discretization (row)
  tmesh = linspace(0, t_final, len_tmesh)'; % Time discretization (column)
  
  % First output argument is max_step_size (ignored). TODO: Move to end of argument
  % list.
  [~, u_true_0, mu_meas, w_meas, g, s_star, s_ini, a_ini] = initial_setup(tmesh, xmesh, use_synthetic_data, initial_data_parameter_s, initial_data_parameter_a);
  
  % Initialize svals and avals
  s_old = s_ini(tmesh);
  a_old = a_ini(tmesh);
  
  % Last argument is functional value (ignored)
  % Second-to-last argument is \tilde{u}(y,t) (ignored)
  [au_xx_S, u_x_S, u_S, u_T, ~, ~] = ...
    Functional(xmesh, tmesh, s_old, a_old, g, u_true_0, s_star, w_meas, mu_meas);

  % Calculate solution of adjoint problem
  % Last argument is psi (ignored)
  [psi_t_S, psi_x_S, psi_S, ~] = ...
    Adjoint(xmesh, tmesh, s_old, a_old, u_T, w_meas, u_S, mu_meas);
  
  % Calculate update direction vector s_update
  s_update = grad_s(tmesh, s_old, w_meas, u_T, mu_meas, u_x_S, psi_x_S, psi_t_S, psi_S, u_S, au_xx_S, s_star);
  if norm(s_update) > norm_update_threshold % Only normalize if the update vector has nonzero norm.
      s_update = s_update/ norm(s_update);
  end
  
  s_update_before = s_update;
  
  % Preconditioning for s(t) gradient
  s_update_after = precond(tmesh, L_s, s_update);
  
  test_grad_s_datafile = 'test_grad_s.mat';
  
  if generate_validation
      save(test_grad_s_datafile, 'tmesh', 's_update_before', 's_update_after', '-mat');
  end
  
  s_update_validation = load(test_grad_s_datafile);
  validation_tmesh = s_update_validation.tmesh;
  validation_s_before = s_update_validation.s_update_before;
  validation_s_after = s_update_validation.s_update_after;
  
  s_update_before_validated = interp1(validation_tmesh, validation_s_before, tmesh);
  s_update_after_validated = interp1(validation_tmesh, validation_s_after, tmesh);
  
  error_out_grad_s = trapz(tmesh, (s_update_before - s_update_before_validated).^2);
  error_out_precond_s = trapz(tmesh, (s_update_after - s_update_after_validated).^2);
end

