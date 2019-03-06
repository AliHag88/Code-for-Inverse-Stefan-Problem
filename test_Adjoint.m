function psiErrorOut = test_Adjoint(len_xmesh, len_tmesh, oscillation)
  % test_Adjoint: Test for correct functioning of Adjoint function.
  % Input Argument:
  %    - len_xmesh: Number of space grid points
  %    - len_tmesh: Number of time grid points
  %    - oscillation: Amount to vary the data for Adjoint inputs.
  %      The expectation for the adjoint problem is that as the boundary data vanishes, the adjoint vanishes.
  if ~exist('oscillation', 'var')
    oscillation = 0;
  end

  xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
  tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

  s_new = @(t) 2 * acoth(exp(exp(t)-1)*coth(1/2));

  u_true = @(x,t) exp(t).*cosh(x);
  u_true_T = @(x) u_true(x, 1);

  boundary_values = s_new(tmesh);
  s_der = est_deriv(boundary_values, tmesh);

  % Set up problem with vanishing boundary data
  u_T = u_true_T(xmesh) + oscillation;
  w_meas = u_true_T;

  mu_meas = @(t) u_true(s_new(t), t);
  u_S = mu_meas(tmesh) + oscillation;

  avals = ones(size(tmesh));

  % Run solver
  [psi_t_S, psi_x_S, psi_S, psi_T, psi] = ...
    Adjoint(xmesh, tmesh, boundary_values, avals, u_T, w_meas, s_der, u_S, mu_meas);

  % Check size of outputs
  assert (all(size(psi) == [len_tmesh, len_xmesh]));

  assert (length(psi_T) == len_xmesh);
  assert (isrow(psi_T));
  assert (length(psi_t_S) == len_tmesh);
  assert (iscolumn(psi_t_S));
  assert (length(psi_x_S) == len_tmesh);
  assert (iscolumn(psi_x_S));
  assert (length(psi_S) == len_tmesh);
  assert (iscolumn(psi_S));

  psiErrorOut = norm(psi, 'fro')*(xmesh(2)-xmesh(1))*(tmesh(2)-tmesh(1));
end
