function visualization(xmesh, tmesh, svals, avals, u, k, J, pausetime)%#ok<INUSL>
  if ~exist('pausetime', 'var')
    pausetime = 0;
  else
    oldState = pause('on');
  end

  len_xmesh = length(xmesh);
  len_tmesh = length(tmesh);

  % Grab analytic data from true_solution
  [~, ~, ~, ~, ~, s_true, ~, ~] = true_solution(tmesh);
  
  % Grab initial data from initial_setup
  [~, ~, ~, ~, ~, ~, s_ini] = initial_setup(tmesh);
  
  %% Surface plot of u on rectangular/square domain
  % subplot(2,3,1)
  % surf(xmesh, tmesh, u)
  % title('Numerical solution u_k (x,t) on square domain.')
  % xlabel('Distance x')
  % ylabel('Time t')

  %% Below, we create an "image" plot, where u-values are translated to different colors.
  subplot(2,3,1)
  imagesc('XData', xmesh, 'YData', tmesh, 'CData', u)
  title('Numerical solution u_k (x,t) on square domain.')
  xlabel('Distance x')
  ylabel('Time t')
  colorbar
  axis tight

  % Create boundary function corresponding to svals
  s_new = @(t) interp1(tmesh, svals, t);

  % u(x,t) visualisation on a non-rectangular domain is based on 2d interpolation.
  [X, T] = meshgrid(xmesh, tmesh);

  % u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
  uf = @(x,t) interp2(X, T, u, x/s_new(t), t, 'linear', NaN);

  % Define new grid on state vector's "native" domain.
  x_new = linspace(0, max(svals), len_xmesh);

  % Initialize and evaluate u_visual, the values on the "native" domain.
  u_visual = zeros(len_tmesh, len_xmesh);
  for i = 1:len_tmesh
    u_visual(i, :) = uf(x_new, tmesh(i));
  end

  %% "image" plot of solution on its "native" domain.
  subplot(2,3,2)
  imagesc('XData', x_new, 'YData', tmesh, 'CData', u_visual)
  hold on
  plot(s_true(tmesh), tmesh, '*', 'color', 'red')
  plot(svals, tmesh, 'O', 'color', 'green')
  plot(s_ini(tmesh), tmesh, '--', 'color', 'white')
  title('True solution u(x,t), initial s(t), true s(t) and s_k(t).')
  xlabel('Distance x')
  ylabel('Time t')
  colorbar
  axis tight

  %% Surface plot of solution on "native" domain.
  % subplot(2,3,4)
  % surf(x_new,tmesh,u_visual)
  % title('Surface graph of u_k (x,t)')
  % xlabel('Distance x')
  % ylabel('Time t')


  %% Create plot of functional values
  subplot(2,3,4)
  scatter(k, J(k), 'filled', 'MarkerFaceColor', 'black')
  title('Cost Functional J')
  xlabel('Iteration k')
  ylabel('Functional Value J(k)')

  if pausetime > 0
    pause(pausetime);
    pause(oldState);
  end
end
