function psiErrorOut = test_Adjoint(len_xmesh, len_tmesh)
% len_xmesh: Number of space grid points
% len_tmesh: Number of time grid points

xmesh = linspace(0,1,len_xmesh); % Space discretization for forward problem
tmesh = linspace(0,1,len_tmesh); % Time discretization for forward problem

s_new = @(t) (3*exp(-t.^2)-1)/2;

boundary_values = s_new(tmesh);
s_der = deriv_est(boundary_values, tmesh);

% Set up problem with vanishing boundary data
u_T = ones(size(xmesh));
w_meas = @(x) ones(size(x));

u_S = -ones(size(tmesh));
mu_meas = @(t) -ones(size(t));

avals = ones(size(tmesh));

% Run solver
[psi_T, psi_t, psi_t_S, psi_S, psi_x_S, psi_x, psi] = ...
    Adjoint(xmesh, tmesh, boundary_values, avals, u_T, w_meas, s_der, u_S, mu_meas);

% Check size of outputs
assert (all(size(psi) == [len_tmesh, len_xmesh]))
assert (all(size(psi_t) == [len_tmesh, len_xmesh]))
assert (all(size(psi_x) == [len_tmesh, len_xmesh]))

assert (length(psi_T) == len_xmesh)
assert (length(psi_t_S) == len_tmesh)
assert (length(psi_x_S) == len_tmesh)
assert (length(psi_S) == len_tmesh)

psiErrorOut = norm(psi, 'fro')*(xmesh(2)-xmesh(1))*(tmesh(2)-tmesh(1));

end
