function visualization(xmesh, tmesh, svals, avals, u, k, J, pausetime, initial_data_parameter_s, initial_data_parameter_a)
  if ~exist('pausetime', 'var')
    pausetime = 0;
  else
    oldState = pause('on');
  end

  len_xmesh = length(xmesh);
  len_tmesh = length(tmesh);

  % Grab analytic data from true_solution
  [~, ~, ~, ~, ~, s_true, u_true, a_true] = true_solution(tmesh);

  % Grab initial data from initial_setup. Set use_synthetic_data=false,
  % since we aren't using the corresponding output parameters.
  [~, ~, ~, ~, ~, ~, s_ini, a_ini] = initial_setup(tmesh, xmesh, false, initial_data_parameter_s, initial_data_parameter_a);

  % Define new grid on state vector's "native" domain.
  x_new = linspace(0, max(svals), len_xmesh);

  % Initialize and evaluate u_visual, the values on the "native" domain.
  u_visual = zeros(len_tmesh, len_xmesh);
  u_true_visual = zeros(len_tmesh, len_xmesh);
  for i = 1:len_tmesh
      sk_curr = svals(i);
      t_curr = tmesh(i);
      s_curr = s_true(t_curr);
      shat_curr = min(sk_curr, s_curr);
      
      u_true_visual(i, :) = u_true(x_new, t_curr);
      u_visual(i, :) = interp1(xmesh * sk_curr, u(i, :), x_new, 'linear', 0);
      
      % Both u_true_visual and u_visual are continued by zero
      u_true_visual(i, x_new > shat_curr) = 0;
      u_visual(i, x_new > shat_curr) = 0;
  end

  subplot(2,3,1);
  %% Image plot of solution on its "native" domain.
  imagesc('XData', x_new, 'YData', tmesh, 'CData', u_visual);
  hold on
  plot(s_true(tmesh), tmesh, '*', 'color', 'red');
  plot(svals, tmesh, 'O', 'color', 'green');
  plot(s_ini(tmesh), tmesh, '--', 'color', 'white');
  title('u_k(x,t), s_{true}(t), s_k(t), and s_{ini}(t).');
  legend({'s_{true}', 's_k', 's_{ini}'}, 'Location', 'best', 'FontSize', 10);
  xlabel('Distance x');
  ylabel('Time t');
  colormap(hot(20)); % Take 20 samples
  colorbar;
  axis([min(x_new), max(x_new), min(tmesh), max(tmesh)]);
  hold off;
  %% End Image plot of solution on "native" domain


  subplot(2,3,2);
  %% Scatter plot of functional values
  plot(J(1:k), '*', 'MarkerFaceColor', 'black');
  title('Cost Functional J');
  xlabel('Iteration k');
  ylabel('Functional Value J(k)');
  axis([0, k, 0, 2*max(J)]);
  %% End plot of functional values


  subplot(2,3,3);
  %% Plot of solution declination
  imagesc('XData', x_new, 'YData', tmesh, 'CData', u_visual - u_true_visual);
  title('u_k - u_{true}');
  xlabel('Distance x');
  ylabel('Time t');
  colormap(hot(20)); % Take 20 samples
  colorbar;
  axis([min(x_new), max(x_new), min(tmesh), max(tmesh)]);
  %% End plot of solution declination in rectangular domain


  subplot(2,3,4)
  %% Line plot of s_true, svals, and s_ini
  plot(tmesh, s_true(tmesh));
  hold on
  plot(tmesh, svals, '*');
  plot(tmesh, s_ini(tmesh), '-');
  legend({'s_{true}', 's_{opt}', 's_{ini}'}, 'Location', 'best', 'FontSize', 10);
  xlabel('Time t');
  ylabel('Position s(t)');
  title('Optimization Output, x=s(t)');
  axis([min(tmesh), max(tmesh), 0, Inf]);
  hold off
  %% End plot of s_true, svals, and s_init


  subplot(2,3,5)
  %% Line plot of a_true, avals, and a_ini
  plot(tmesh, a_true(tmesh));
  hold on % Must go after first plot on axes
  plot(tmesh, avals, '*');
  plot(tmesh, a_ini(tmesh), '-');
  legend({'a_{true}', 'a_{opt}', 'a_{ini}'}, 'Location', 'best', 'FontSize', 10);
  xlabel('Time t');
  ylabel('Value a(t)');
  title('Optimization Output, a(t)');
  axis([min(tmesh), max(tmesh), 0, Inf]);
  hold off
  %% End plot of a_true, avals, and a_init


  if pausetime > 0
    pause(pausetime);
    pause(oldState);
  end
end

%% Spare plotting "recipes"

%% Surface plot of u on rectangular/square domain
% surf(xmesh, tmesh, u)
% title('Numerical solution u_k (x,t) on square domain.')
% xlabel('Distance x')
% ylabel('Time t')
%% End surface plot of u on rectangular domain

%% Surface plot of solution on "native" domain.
% surf(x_new,tmesh,u_visual)
% title('Surface graph of u_k (x,t)')
% xlabel('Distance x')
% ylabel('Time t')
%% End surface plot of solution on "native" domain

%% "image" plot, where u-values are translated to different colors.
%imagesc('XData', xmesh, 'YData', tmesh, 'CData', u)
%title('Numerical solution u_k (x,t) on square domain.')
%xlabel('Distance x')
%ylabel('Time t')
%colorbar
%axis tight
%% End Image plot of solution on rectangular domain
