xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

[ ...
  u_true_T, u_true_0, u_true_S, ...
  aux_true_0, ...
  s_star, s_true, ...
  u_true, a_true ...
] = true_solution(tmesh);

u_true_fullDomain = @(x,t) u_true(x,t).*(x<s_true(t))

% Define grid for Forward problem solution
[X, T] = meshgrid(xmesh, tmesh);

% Easiest to index U, X, and T with a common linear index (u_true_fullDomain
% is not currently vectorized)
for i in 1:numel(X)
  U[i] = u_true_fullDomain(X(i), T(i));
end

subplot(1,2,1);
surf(X, T, U);
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

