function psiErrorOut = test_Adjoint(len_xmesh, len_tmesh, do_visualization, use_idealized_data, oscillation)
  % test_Adjoint: Test for correct functioning of Adjoint problem solver
  % Input Arguments:
  %    - len_xmesh: Number of space grid points
  %    - len_tmesh: Number of time grid points
  
  if ~exist('do_visualization', 'var')
      do_visualization = false;
  end
  if ~exist('use_idealized_data', 'var')
      use_idealized_data = true;
  end
  if ~exist('oscillation', 'var')
      oscillation = 0;
  end

  xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
  tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem
  
  [u_true_T, u_true_0, u_true_S, g, ~, s_true, ~, a_true] = true_solution(tmesh);

  boundary_values = s_true(tmesh);

  avals = a_true(tmesh);
  mu_meas = u_true_S;
  w_meas = u_true_T;
  
  if use_idealized_data
      u_S = u_true_S(tmesh) + oscillation;
      u_T = @(x) u_true_T(x) + oscillation;
  else
      % Set up problem with analytic solution as input for adjoint.
      [~, ~, u_S, u_T, ~] = Forward(xmesh, tmesh, boundary_values, avals, g, u_true_0);
  end
  
  % Run adjoint solver
  [psi_t_S, psi_x_S, psi_S, psi] = ...
    Adjoint(xmesh, tmesh, boundary_values, avals, u_T, w_meas, u_S, mu_meas);

  % Check size of outputs
  assert (all(size(psi) == [len_tmesh, len_xmesh]));
  assert (length(psi_t_S) == len_tmesh);
  assert (iscolumn(psi_t_S));
  assert (length(psi_x_S) == len_tmesh);
  assert (iscolumn(psi_x_S));
  assert (length(psi_S) == len_tmesh);
  assert (iscolumn(psi_S));

  psiErrorOut = norm(psi, 'fro')*(xmesh(2)-xmesh(1))*(tmesh(2)-tmesh(1));
  
  if do_visualization
      figure();
      subplot(1,4,1)
      xmesh_final = xmesh * boundary_values(end);
      plot(xmesh_final, psi(end, :));
      xlabel('x');
      ylabel('\psi');
      title('Value of \psi at t=T');
      
      subplot(1,4,2);
      xmesh_initial = xmesh * boundary_values(1);
      plot(xmesh_initial, psi(1, :));
      xlabel('x');
      ylabel('\psi');
      title('Value of \psi at t=0');
      
      subplot(1,4,3);
      plot(tmesh, psi_S);
      xlabel('t');
      ylabel('\psi');
      title('Value of \psi at x=s(t)');
      
      subplot(1,4,4);
      imagesc('XData', xmesh, 'YData', tmesh, 'CData', psi);
      xlabel('x');
      ylabel('t');
      title('\psi');
      colormap(hot(20)); % Take 20 samples
      colorbar();
      axis([min(xmesh), max(xmesh), min(tmesh), max(tmesh)])
  end
end
