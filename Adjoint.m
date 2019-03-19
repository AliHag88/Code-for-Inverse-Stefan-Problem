function [psi_t_S, psi_x_S, psi_S, psi_T, psi] = Adjoint(xmesh,tmesh,svals,avals,u_T,w_meas,s_der,u_S,mu_meas)
  % Adjoint: Compute values of adjoint for ISP
  % Input Arguments:
  %    - xmesh: Space discretization
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points
  %    - u_T: Vector of values u(x,T;v) after transformation y=x/s(t)
  %    - w_meas: Function of x, representing measurement u(x,T)=:w(x)
  %    - s_der: Vector of values s'(t)
  %    - u_S: Vector of trace values u(s(t),t;v)
  %    - mu_meas: Function of t, representing measurement u(s(t),t)=:mu(t)

%%%
%%% Set up solver/parameters
%%%

  t_final = tmesh(end);

  % Interpolate s(t) values to form \bar{s}(t)=s(T-t)
  s_bar = @(t) interp1(tmesh, svals, t_final - t);

  % Interpolate s' values to form s'(T-t)
  sder = @(t) interp1(tmesh, s_der, t_final - t);

  % Interpolate a values to form \bar{a}(t) = a(T-t)
  af = @(t) interp1(tmesh, avals, t_final - t);

  % u_T is not transformed to a non-rectangular domain by Forward.m
  uT = @(y) interp1(xmesh, u_T, y);

  uS = @(t) interp1(tmesh, u_S, t_final - t);

  mu = @(t) mu_meas(t_final - t);

%%%
%%% Calculate solution on rectangular domain. These represent \tilde{\psi}(y,t)=psi(ys(t),t)
%%%
psi = pdeSolver(xmesh, tmesh, af, s_bar, sder, uT, w_meas, uS, mu);

%%%
%%% Calculate values derived from psi
%%%

% spatial stepsize
h = xmesh(2) - xmesh(1);

% time stepsize
tau = tmesh(2) - tmesh(1);

% Values of psi(y,t) at y=1, which is psi(s(t),t).
psi_S = psi(:, end);

% Values of psi(x,t) at t=T. Note that these values are not transformed to
% the non-rectangular domain.
psi_T = psi(end, :);

% In order to compute psi_x(x,t) from the output from pdeSolver, which
% corresponds to \tilde{psi}(y,t)=psi(ys(t),t), we need to return
% psi_x(s(t),t) = \tilde{\psi}_y(1,t)/s(t) = psi_x(:, end) ./ svals
psi_x = est_x_partial(psi, h);
psi_x_S = psi_x(:, end) ./ svals;

% TODO: Compute these using second derivative information through PDE.
% This would allow Adjoint and Forward to use nearly the same post-processing
% step.
% In a similar way to above, we would like to produce values of psi_t(x,t),
% but pdeSolver produces values \tilde{psi}(y,t)=psi(ys(t),t)
% Hence we return
% psi_t(x,t) = - x s'(t) / s^2(t) \tilde{\psi}_y(y,t) + \tilde{\psi}_t(y,t)
%            = - y s'(t)/s(t) \tilde{\psi}_y(y,t) + \tilde{\psi}_t(y,t)
% We first compute \tilde{psi}_t(y,t):
psi_t = est_t_partial(psi, tau);
% We need only the trace at x=s(t), which corresponds to y=1
psi_t_S = -s_der .* psi_x(:, end) + psi_t(:, end);

end
% End Main Function


% Begin Subfunctions
function psi = pdeSolver(xmesh, tmesh, af, sf, sder, uT, uTrueT, uS, mU)
  m = 0; % Symmetry of the problem. 0=slab (rectangular)

  % PDE to be solved is
  %    c u_t = x^{-m} (x^m f(x,t,u, DuDx))_x + s
  % The output values below give the corresponding value at each point (x, t, u, u_x)
  pde = @(x,t,u, DuDx) ...
         deal(...
          sf(t)^2, ... #c
          af(t)*DuDx, ... #f
          sf(t)*sder(t)*x*DuDx ... #s
         );

  % Initial condition u(x,t_0)
  ic = @(x) 2*(uT(x)-uTrueT(x));

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
              0, ... # pl
              1, ... # ql
              sf(t)*(-sder(t)*ur+2*(uS(t)-mU(t))), ... # pr
              1 ... # qr
            );
  
  % Calculate solution evaluated at each point of xmesh and tmesh.
  sol = pdepe(m, pde, ic, bc, xmesh, tmesh);

  % Extract the first solution component as psi.
  % Note: This solution is defined on tmesh * xmesh.
  % Note: This last function reverses time so psi is "defined" on (0,1)x(0,T)
  psi = flipud(sol(:, :, 1));
end
% End Subfunctions
