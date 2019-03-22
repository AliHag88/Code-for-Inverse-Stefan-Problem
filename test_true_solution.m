function [tf] = test_true_solution(len_xmesh, len_tmesh, visualize)
if ~exist('len_xmesh', 'var')
    len_xmesh = 100;
end
if ~exist('len_tmesh', 'var')
    len_tmesh = 100;
end
if ~exist('visualize', 'var')
    visualize = false;
end

tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

[ ...
  u_true_T, u_true_0, u_true_S, ...
  aux_true_S, ...
  s_star, ...
  s_true, u_true, ...
  a_true ...
  ] = true_solution(tmesh);

boundary_values = s_true(tmesh);
xmesh = linspace(0, max(boundary_values), len_xmesh);

% Define grid for Forward problem solution
[X, T] = meshgrid(xmesh, tmesh);

% Easiest to index U, X, and T with a common linear index (u_true_fullDomain
% is not currently vectorized)
U = zeros(size(X));
for i = 1:numel(X)
    if X(i) > s_true(T(i))
        continue
    end
    U(i) = u_true(X(i), T(i));
end

% Check for correct shape and size of output arguments (without more care
% this _can_ go wrong!
assert(all(U(:) <= 0))

assert(isrow(u_true_T(xmesh)))
assert(numel(u_true_T(xmesh)) == len_xmesh)

assert(isrow(u_true_0(xmesh)))
assert(numel(u_true_0(xmesh)) == len_xmesh)

assert(iscolumn(u_true_S(tmesh)))
assert(numel(u_true_S(tmesh)) == len_tmesh)

assert(iscolumn(aux_true_S(tmesh)))
assert(numel(aux_true_S(tmesh)) == len_tmesh)

assert(iscolumn(s_true(tmesh)))
assert(numel(s_true(tmesh)) == len_tmesh)

assert(iscolumn(a_true(tmesh)))
assert(numel(a_true(tmesh)) == len_tmesh)

assert(isscalar(s_star))
assert(s_star > 0)

if visualize
  % For graphing, remove values if they are close to zero or positive
  U(U>-1e-8) = NaN;

  figure();

  subplot(1,2,1);
  surf(X, T, U, 'EdgeColor', 'none');
  xlabel('x');
  ylabel('t');
  zlabel('u');
  title('Surface plot of Analytical Solution');

  subplot(1,2,2);
  hold on
  imagesc('XData', xmesh, 'YData', tmesh, 'CData', U);
  plot(s_true(tmesh), tmesh, 'color', 'red');
  xlabel('Distance x');
  ylabel('Time t');
  title('Analytical Solution and Boundary Curve');
  colorbar
  axis tight
  hold off

  figure()
  plot(tmesh, u_true_S(tmesh))
  title('u_{true}(s(t),t) == 0')
end

tf = true;
end