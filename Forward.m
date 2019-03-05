function [au_xx_S, u_x_S, u_S, u_T, u] = Forward(xmesh, tmesh, svals, avals, uInitial)
  % Forward: Solve forward problem on given mesh with given data.
  % Input Arguments:
  %    - xmesh: Space discretization on [0,1]
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points
  %    - uInitial: Function u(x,0)
  %
  % Output Arguments:
  %    - au_xx_S: Vector of trace values (a u_x)_x on x=s(t)
  %    - u_x_S: Vector of trace values u_x on x=s(t)
  %    - u_S: Vector of trace values u on x=s(t)
  %    - u_T: Vector of trace values u on t=T
  %    - u: Matrix of output values. Rows are time-levels, columns are space-levels

%%%
%%% Set up solver/parameters
%%%

% Define the vector of derivatives for s values
s_der = deriv_est(svals, tmesh);

% Interpolate s values to form s(t)
s_ini = @(t) interp1(tmesh, svals, t);

% Interpolate a values to form a(t)
af = @(t) interp1(tmesh, avals, t);

% Interpolate s' values to form s'(t)
sder = @(t) interp1(tmesh, s_der, t);

%%%
%%% Calculate solution on rectangular domain
%%%
u = pdeSolver(xmesh, tmesh, af, s_ini, sder, uInitial);

%%%
%%% Transform solution to non-rectangular domain
%%%
% Define "new" (non-rectangular) boundary curve
% s_new =@(t) interp1(tmesh, svals, t);

% u(x,t) visualisation on a non-rectangular domain is based on 2d interpolation.
% [X, T] = meshgrid(xmesh, tmesh);

% u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
% uf = @(x,t) interp2(X, T, u, x/s_new(t), t, 'linear', NaN);

%%%
%%% Calculate values of solution at requested locations
%%%

%% Values Of u(x,t) at t=T
u_T = u(end, :);

%% Values of u(x,t) at s(t)
u_S = u(:, end);

% Number of space grid points
len_xmesh = length(xmesh);

% spatial stepsize
h = xmesh(3)-xmesh(2);

%% u_x (x,t) matrix of derivatives of u(x,t)
u_x = zeros(size(u));

% Use forward difference for internal points
u_x(:, 1:end-1) = u(:, 2:end) - u(:, 1:end-1);

% Use backward difference at the right-hand side
u_x(:, end) = u(:, end) - u(:, end - 1);

% Scale all at once
u_x = u_x / h;

for i = 1:len_xmesh 
  u_x(:,i) = u_x(:,i) ./ svals(:);
end


% Values of u_x (x,t) at s(t)
u_x_S = u_x(:, end);

%% Values of a(t)u_x(x,t)
au_x = zeros(size(u));
for i = 1:len_xmesh
  au_x(:, i) = avals(:) .* u_x(:, i);
end


%% values of (a(t)u_x(x,t))_x
au_xx = zeros(size(u));

% Calculate values inside domain using forward difference
au_xx(:, 1:end-1) = au_x(:, 2:end) - au_x(:, 1:end-1);

% Calculate values on right-hand side of domain using backward difference
au_xx(:, end) = au_xx(:, end-1);

% Scale all differences at once
au_xx = au_xx / h;

for i = 1:len_xmesh
  au_xx(:,i) = au_xx(:,i) ./ ((svals(:)).^2);
end

% values of (a(t)u_x(x,t))_x at s(t)
au_xx_S = au_xx(:, end);

end
%%% End Main Function

%%% Begin Subfunctions
function u = pdeSolver(xmesh, tmesh, af, sf, sder, uTrue0)
    m = 0; % Symmetry of the problem. 0=slab (rectangular)

    % PDE to be solved is
    %    c u_t = x^{-m} (x^m f(x,t,u, DuDx))_x + s
    % The output values below give the corresponding value at each point (x, t, u, u_x)
    pde = @(x,t,u,DuDx) ...
           deal(...
                 sf(t).^2, ... # c
                 af(t)*DuDx, ... # f
                 sf(t)*sder(t)*x*DuDx ... # s
               );

    % Initial condition u(x,t_0)
    ic = uTrue0;

    % Boundary condition at x=xl and x=xr.
    % Suffix l corresponds to left-hand side of domain,
    % Suffix r corresponds to right-hand side of domain.
    % Boundary condition takes the form
    %    p(x,t,u) + q(x,t) f(x,t,u,u_x) = 0
    % where f is defined below.
    % In particular, the inputs correspond to (x,u) at x=xl and (x,u) at x=xr
    % and the time t, and the outputs correspond to (p,q) at x=xl and (p,q) at
    % x=xr
    bc = @(xl, ul, xr, ur, t) ...
          deal(...
                -sf(t), ... # pl
                1, ... # ql
                sder(t)*sf(t), ... # pr
                1 ... # qr
              );

    % Calculate solution evaluated at each point of xmesh and tmesh.
    sol = pdepe(m, pde, ic, bc, xmesh, tmesh);

    % Extract the first solution component as u.
    % Note: This solution is defined on tmesh * xmesh.
    % TODO: Transform to function u(x,t) within this routine.
    u = sol(:, :, 1);
end
%%% End Subfunctions
