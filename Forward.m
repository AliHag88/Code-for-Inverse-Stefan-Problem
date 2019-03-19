function [au_xx_S, u_x_S, u_S, u_T, u] = Forward(xmesh, tmesh, svals, avals, g, uInitial)
  % Forward: Solve forward problem on given mesh with given data.
  % Input Arguments:
  %    - xmesh: Space discretization on [0,1]
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points
  %    - g: Function a(t) u_x(0, t) =: g(t)
  %    - uInitial: Function u(x,0)
  %
  % Output Arguments:
  %    - au_xx_S: Vector of trace values (a u_x)_x on x=s(t)
  %    - u_x_S: Vector of trace values u_x on x=s(t)
  %    - u_S: Vector of trace values u on x=s(t)
  %    - u_T: Vector of trace values u on t=T after transformation x/s(t)
  %    - u: Matrix of output values. Rows are time-levels, columns are space-levels

%%%
%%% Set up solver/parameters
%%%

% Define the vector of derivatives for s values
s_der = est_deriv(svals, tmesh);

% Interpolate s values to form s(t)
s_current = @(t) interp1(tmesh, svals, t);

% Interpolate a values to form a(t)
af = @(t) interp1(tmesh, avals, t);

% Interpolate s' values to form s'(t)
sder = @(t) interp1(tmesh, s_der, t);

uInitial_tilde = @(y) uInitial(y*svals(1));

%%%
%%% Calculate solution u(y,t) on rectangular domain
%%%
u = pdeSolver(xmesh, tmesh, af, s_current, sder, g, uInitial_tilde);

%%%
%%% Calculate values of solution at requested locations
%%%

%% Values of u(y,t) at t=T. Note that these values are not transformed to
% the non-rectangular domain.
u_T = u(end, :);

%% Values of u(y,t) at y=1, which is u(s(t), t).
u_S = u(:, end);

% spatial stepsize
h = xmesh(2) - xmesh(1);

% In order to compute u_x(x,t) from the output from pdeSolver, which
% corresponds to \tilde{u}(y,t)=psi(ys(t),t), we need to return
% psi_x(s(t),t) = \tilde{u}_y(1,t)/s(t) = u_x(:, end) ./ svals
u_x = est_x_partial(u, h);
u_x_S = u_x(:, end) ./ svals;

% Similarly, we compute a(t) u_{xx}(s(t), t) = a(t)\tilde{u}_{yy}(1,t)/s^2(t)
u_xx = est_x_partial(u_x, h);
au_xx_S = (avals .* u_xx(:, end)) ./ (svals.^2);

end
%%% End Main Function

%%% Begin Subfunctions
function u = pdeSolver(xmesh, tmesh, af, sf, sder, g, uTrue0)
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
                g(t)*sf(t), ... # pl
                1, ... # ql
                -sder(t)*sf(t), ... # pr
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
