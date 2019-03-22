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
  %    - u_T: Function u(x, T)

%%%
%%% Set up solver/parameters
%%%

% Define the vector of derivatives for s values, and interpolate s' values
% to form s'(t)
sder = @(t) interp1(tmesh, est_deriv(svals, tmesh), t);

%%%
%%% Calculate solution u(y,t) on rectangular domain
%%%
u = pdeSolver(...
    xmesh, tmesh, svals, avals, sder, ...
    g, ... % aux_0
    @(t) zeros(size(t)), ... robin_coeff
    sder, ... robin_rhs
    uInitial);

%%%
%%% Calculate values of solution at requested locations
%%%

%% Derive trace at final moment
s_curr = svals(end);
xmesh_curr = xmesh * s_curr;
[u_T_vals, ~] = pdeval(0, xmesh_curr, u(end, :), xmesh_curr); % second output is u_x(x, T)
u_T = @(x) interp1(xmesh_curr, u_T_vals, x/s_curr, 'linear', 'extrap');

%% Calculate traces of solution using pdeval
u_S = zeros(size(tmesh));
u_x_S = zeros(size(tmesh));
au_xx_S = zeros(size(tmesh));

for k = 1:length(tmesh)
    s_curr = svals(k);
    xmesh_curr = xmesh * s_curr;
    [u_Sk, u_x_Sk] = pdeval(0, xmesh_curr, u(k, :), [xmesh_curr(end - 1), xmesh_curr(end)]);
    u_S(k) = u_Sk(2);
    u_x_S(k) = u_x_Sk(2) / s_curr;
    au_xx_S(k) = diff(u_x_Sk) / ((xmesh_curr(end) - xmesh_curr(end - 1))*s_curr^2);
end
end
%%% End Main Function

% Alternative implementation using BVP solver
% function u = pdeSolver(xmesh, tmesh, af, sf, sder, g, uTrue0)
%     tau = tmesh(2) - tmesh(1);
%     u = zeros(length(tmesh), length(xmesh));
%     
%     % Save initial data
%     k = 1;
%     t_curr = tmesh(k);
%     s_curr = sf(t_curr);
%     xmesh_curr = xmesh * s_curr;
% 
%     u(k, :) = arrayfun(uTrue0, xmesh_curr);
%     
%     % Create interpolation of initial data
%     uPrev = @(x) interp1(xmesh_curr, u(k, :), x, 'linear', 'extrap');
%     for k = 2:length(tmesh)
%         t_curr = tmesh(k);
%         a_curr = af(t_curr);
%         s_curr = sf(t_curr);
%         g_curr = g(t_curr);
%         xmesh_curr = xmesh * s_curr;
%         
%         rhs = @(x, y) [
%             y(2);
%             (y(1) - uPrev(x))/(tau * a_curr);
%             ];
%         bcs = @(ya, yb) [
%             a_curr*ya(2) - g_curr;
%             a_curr*yb(2) - sder(t_curr)
%             ];
%         guess = @(x)[
%             uPrev(x);
%             -g_curr/a_curr
%             ];
%         solinit = bvpinit(xmesh_curr, guess);
%         sol = bvp4c(rhs, bcs, solinit);
%                 
%         uPrev = @(x) interp1(sol.x, sol.y(1, :), x, 'linear', 'extrap');
%         u(k, :) = arrayfun(uPrev, xmesh_curr);
%     end
% end
%%% End Subfunctions
