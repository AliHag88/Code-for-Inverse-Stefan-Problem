function u = pdeSolver(xmesh, tmesh, svals, avals, sder, aux_0, robin_coeff, robin_rhs, u_initial)
  % pdeSolver: Solve PDE problem non-rectangular domain with given data.
  % Solves the problem
  % (1)    u_t = (a(t) u_x)_x , 0<x<s(t), 0<t<T
  % (2)    u(x,0) = u_initial(x), 0<x<s(0)
  % (3)    a(t)u_x(0,t) = aux_0(t)
  % (4)    a(t)u_x(s(t),t) + r(t) u(s(t),t) = p(t)
  %
  % on a given time grid by transforming the domain to [0,1]x[0,T].
  %
  % Input Arguments:
  %    - xmesh: Space discretization on [0,1]
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points from
  %             condition (1).
  %    - aux_0: Function aux_0(t) in condition (3)
  %    - robin_coeff: Function r(t) in condition (4)
  %    - robin_rhs: Function p(t) in condition (4)
  %    - u_initial: Function u(x,0) in condition (2)
  %
  % Output Arguments:
  %    - u: Matrix with length(tmesh) rows and length(xmesh) columns.
  % On each time level u(i, :), the values correspond to the spatial grid
  % xmesh * svals(i).

%%%
%%% Set up solver/parameters
%%%

% Interpolate s values to form s(t)
  sf = @(t) interp1(tmesh, svals, t);

% Interpolate a values to form a(t)
  af = @(t) interp1(tmesh, avals, t);


  m = 0; % Symmetry of the problem. 0=slab (rectangular)

    % PDE to be solved is
    %    c u_t = y^{-m} (y^m f(y,t,u, DuDy))_y + s
    % The output values below give the corresponding value at each point (x, t, u, u_x)
    pde = @(y,t,u,DuDx) ...
           deal(...
                 sf(t), ... # c
                 af(t)*DuDx / sf(t), ... # f
                 sder(t)*y*DuDx ... # s
               );

    % Initial condition u(x,t_0)
    ic = @(y) u_initial(y*svals(1));

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
                -aux_0(t), ... # pl
                1, ... # ql
                robin_coeff(t)*ur-robin_rhs(t), ... # pr
                1 ... # qr
              );

    % Calculate solution evaluated at each point of xmesh and tmesh.
    sol = pdepe(m, pde, ic, bc, xmesh, tmesh);

    % Extract the first solution component as u.
    % Note: This solution is defined on tmesh * xmesh.
    u = sol(:, :, 1);
end
