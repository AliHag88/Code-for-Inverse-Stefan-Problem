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

  [~, ~, mu_meas, w_meas, ~, ~] = initial_setup(tmesh, xmesh);
  
  [~, u_true_0, ~, g, ~, s_true, ~, a_true] = true_solution(tmesh);

  boundary_values = s_true(tmesh);
  s_der = est_deriv(boundary_values, tmesh);

  avals = a_true(tmesh);

  % Set up problem with analytic solution as input for adjoint.
  [~, ~, u_S, u_T, ~] = Forward(xmesh, tmesh, boundary_values, avals, g, u_true_0);

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
