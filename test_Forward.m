function errorOut = test_Forward(len_xmesh, len_tmesh)
% len_xmesh: Number of space grid points
% len_tmesh: Number of time grid points

xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
tmesh = linspace(0,1,len_tmesh)'; % Time discretization for forward problem

s_new = @(t) 2 * acoth(exp(exp(t)-1)*coth(1/2));

u_true = @(x,t) exp(t).*cosh(x);
u_true_0 = @(x) u_true(x, 0);

boundary_values = s_new(tmesh);
avals = ones(size(tmesh));
g = @(t) zeros(size(t));

[au_xx_S, u_x_S, u_S, u_T, u] = ...
  Forward(xmesh, tmesh, boundary_values, avals, g, u_true_0);

% Define grid for Forward problem solution
[X, T] = meshgrid(xmesh, tmesh);

% u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
u_interp = @(x,t) interp2(X, T, u, x/s_new(t), t, 'linear', NaN);

x_new = linspace(0, max(boundary_values), len_xmesh);

diff_visual = zeros(len_tmesh, len_xmesh);
for i = 1:len_tmesh
    for j = 1:len_xmesh
        if x_new(j) > boundary_values(i)
            break
        end
        diff_visual(i, j) = u_true(x_new(j), tmesh(i)) - u_interp(x_new(j), tmesh(i));
    end
end

assert (all(size(u) == [len_tmesh, len_xmesh]));
assert (length(au_xx_S) == len_tmesh);
assert (iscolumn(au_xx_S));
assert (length(u_S) == len_tmesh);
assert (iscolumn(u_S));
assert (length(u_x_S) == len_tmesh);
assert (iscolumn(u_x_S));
assert (length(u_T) == len_xmesh);
assert (isrow(u_T));

errorOut = norm(diff_visual, 'fro')*(x_new(2)-x_new(1))*(tmesh(2)-tmesh(1));

end
