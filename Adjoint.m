function [psi_T,psi_t,psi_t_S,psi_S,psi_x_S,psi_x,psi] = Adjoint(xmesh,tmesh,svals,avals,u_T,w_meas,s_der,u_S,mu_meas)
  % Adjoint: Compute values of adjoint for ISP
  % Input Arguments:
  %    - xmesh: Space discretization
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points
  %    - u_T: Vector of values u(x,T;v)
  %    - w_meas: Function of x, representing measurement u(x,T)=:w(x)
  %    - s_der: Vector of values s'(t)
  %    - u_S: Vector of trace values u(s(t),t;v)
  %    - mu_meas: Function of t, representing measurement u(s(t),t)=:mu(t)

%%%
%%% Set up solver/parameters
%%%

  t_final = tmesh(end);

  % Interpolate s values to form s(t)
  s_ini = @(t) interp1(tmesh, fliplr(svals), t);

  % Interpolate s' values to form s'(t)
  sder = @(t) interp1(tmesh, fliplr(s_der), t);

  % Interpolate a values to form a(t)
  af = @(t) interp1(tmesh, fliplr(avals), t);

  uT = @(x) interp1(xmesh, u_T, x);

  uS = @(t) interp1(tmesh, fliplr(u_S), t);

  mu = @(t) mu_meas(t_final - t);

%%%
%%% Calculate solution on rectangular domain
%%%
psi = pdeSolver(xmesh, tmesh, af, s_ini, sder, uT, w_meas, uS, mu);

%%%
%%% Calculate values derived from psi
%%%

% Number of space grid points
len_xmesh = length(xmesh);

% spatial stepsize
h = xmesh(3)-xmesh(2);

% time stepsize
tau = tmesh(2)-tmesh(1);

% psi_x (x,t) matrix of derivatives of psi(x,t)
psi_x = zeros(size(psi));
psi_x(:, 1:end-1) = psi(:, 2:end) - psi(:, 1:end-1);
psi_x(:, end) = psi(:, end) - psi(:, end-1);
% Scale all values
psi_x = psi_x / h;

for i = 1:len_xmesh
  psi_x(:,i) = psi_x(:,i) ./ svals(:);
end

% Values of psi_x (x,t) at s(t)
psi_x_S = psi_x(:, end);


% psi_t (x,t) matrix of derivatives of psi(x,t) w.r.t t
psi_t=zeros(size(psi));
psi_t(1:end-1, :) = psi(2:end,:) - psi(1:end-1,:);
psi_t(end, :) = psi(end, :) - psi(end-1,:);
% Scale all values
psi_t = psi_t / tau;

for i = 1:len_xmesh
    psi_t(:,i) = psi_t(:,i) - psi_x(:,i) .* s_der(:);
end

% Values of psi_t (x,t) at s(t)
psi_t_S = psi_t(:, end);

% Values of psi (x,t) at s(t)
psi_S = psi(:, end);

% Values Of psi(x,t) at t=T
psi_T = psi(end, :);
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
          -sf(t)^2, ... #c
          -af(t)*DuDx, ... #f
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
              sder(t)*sf(t)*ul+2*sf(t)*(uS(t)-mU(t)), ... # pr
              -1 ... # qr
            );

  % Calculate solution evaluated at each point of xmesh and tmesh.
  sol = pdepe(m, pde, ic, bc, xmesh, tmesh);

  % Extract the first solution component as psi.
  % Note: This solution is defined on tmesh * xmesh.
  % TODO: Transform to function psi(x,t) within this routine.
  psi = flipud(sol(:, :, 1));
end
% End Subfunctions
