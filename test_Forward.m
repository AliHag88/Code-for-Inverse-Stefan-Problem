function [uErrorOut, u_SErrorOut, u_TErrorOut] = test_Forward(len_xmesh, len_tmesh, do_visualization)
  % test_Forward: Test for correct functioning of Forward problem solver
  % len_xmesh: Number of space grid points
  % len_tmesh: Number of time grid points
  if ~exist('do_visualization', 'var')
      do_visualization = false;
  end

  xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
  tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

  [u_true_T, u_true_0, u_true_S, g, ~, s_true, u_true, a_true] = true_solution(tmesh);

  boundary_values = s_true(tmesh);
  avals = a_true(tmesh);

  [au_xx_S, u_x_S, u_S, u_T, u] = ...
    Forward(xmesh, tmesh, boundary_values, avals, g, u_true_0);

  diff_visual = zeros(len_tmesh, len_xmesh);
  for i = 1:len_tmesh
      s_curr = boundary_values(i);
      for j = 1:len_xmesh
          diff_visual(i, j) = u_true(xmesh(j)*s_curr, tmesh(i)) - u(i, j);
      end
  end
  
  % Initial data should match
  assert(norm(diff_visual(1, :))^2*(xmesh(2)-xmesh(1)) < 1e-3);
  
  % Check size of outputs
  assert (all(size(u) == [len_tmesh, len_xmesh]));
  assert (length(au_xx_S) == len_tmesh);
  assert (iscolumn(au_xx_S));
  assert (length(u_S) == len_tmesh);
  assert (iscolumn(u_S));
  assert (length(u_x_S) == len_tmesh);
  assert (iscolumn(u_x_S));
  assert (length(u_T(xmesh)) == len_xmesh);
  assert (isrow(u_T(xmesh)));

  uErrorOut = norm(diff_visual, 'fro')*(xmesh(2)-xmesh(1))*(tmesh(2)-tmesh(1));
  
  s_final = boundary_values(end);
  xmesh_final = xmesh * s_final;
  u_TErrorOut = norm(u_T(xmesh_final) - u_true_T(xmesh_final))^2 * (xmesh(2) - xmesh(1));
  u_SErrorOut = norm(u_S - u_true_S(tmesh))^2 * (tmesh(2) - tmesh(1));
  
  if do_visualization
      figure();
      subplot(1,3,1);
      plot(xmesh_final, u_T(xmesh_final), xmesh_final, u_true_T(xmesh_final));
      legend({'u(x,T)', 'u_{true}(x,T)'});
      xlabel('x');
      ylabel('u');
      title('Measurements at t=T');
      
      subplot(1,3,2);
      plot(tmesh, u_S, tmesh, u_true_S(tmesh));
      legend({'u(s(t),t)', 'u_{true}(s(t),t)'});
      xlabel('t');
      ylabel('u');
      title('Measurements at x=s(t)');
      
      subplot(1,3,3);
      imagesc('XData', xmesh, 'YData', tmesh, 'CData', diff_visual);
      xlabel('x');
      ylabel('t');
      title('u_{true}-u');
      colormap(hot(20)); % Take 20 samples
      colorbar();
      axis([min(xmesh), max(xmesh), min(tmesh), max(tmesh)])
  end
end
