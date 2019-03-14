function [] = viz_true_solution(len_xmesh, len_tmesh)
if ~exist('len_xmesh', 'var')
    len_xmesh = 100;
end
if ~exist('len_tmesh', 'var')
    len_tmesh = 100;
end

tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

[ ...
  ~, ~, ~, ...
  ~, ...
  ~, s_true, ...
  u_true, ~...
] = true_solution(tmesh);

boundary_values = s_true(tmesh);
xmesh = linspace(0, max(boundary_values), len_xmesh);

u_true_fullDomain = @(x,t) u_true(x,t).*(x<s_true(t));

% Define grid for Forward problem solution
[X, T] = meshgrid(xmesh, tmesh);

% Easiest to index U, X, and T with a common linear index (u_true_fullDomain
% is not currently vectorized)
U = zeros(size(X));
for i = 1:numel(X)
  U(i) = u_true_fullDomain(X(i), T(i));
end
% If values are close to zero, remove them entirely
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

end