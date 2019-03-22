function [psi_t_S, psi_x_S, psi_S, psi] = Adjoint(xmesh,tmesh,svals,avals,u_T,w_meas,u_S,mu_meas)
  % Adjoint: Compute values of adjoint for ISP
  % Input Arguments:
  %    - xmesh: Space discretization
  %    - tmesh: Time discretization
  %    - svals: Position of boundary at time grid points
  %    - avals: Diffusion coefficient a(t) at time grid points
  %    - u_T: Function u(x,T;v)
  %    - w_meas: Function of x, representing measurement u(x,T)=:w(x)
  %    - u_S: Vector of trace values u(s(t),t)
  %    - mu_meas: Function of t, representing measurement u(s(t),t)=:mu(t)

%%%
%%% Set up solver/parameters
%%%

  t_final = tmesh(end);
  
  % Calculate estimate for s' at grid points
  s_der = est_deriv(svals, tmesh);

  % Interpolate s' values to form s'(T-t)
  sder_ = @(t) interp1(tmesh, s_der, t_final - t);

  uS = @(t) interp1(tmesh, u_S, t);

%%%
%%% Calculate solution on rectangular domain. Output represents \tilde{\psi}(y,t)=psi(ys(t),t)
%%%
psi = flipud(pdeSolver(...
    xmesh, tmesh, flipud(avals), flipud(svals), sder_,...
    @(t) zeros(size(t)), ... % aux_0
    @(t) -sder_(t), ... % robin_coeff
    @(t) uS(t_final - t) - mu_meas(t_final - t_final), ... % robin_rhs
    @(x) u_T(x) - w_meas(x) ... % u_initial
    ));

%%%
%%% Calculate values derived from psi
%%%
%% Calculate traces of solution using pdeval
tau = tmesh(2) - tmesh(1);
psi_S = zeros(size(tmesh));
psi_x_S = zeros(size(tmesh));
psi_t_S = zeros(size(tmesh));

for k = 1:length(tmesh)
    s_curr = svals(k);
    [psi_S(k), psi_x_Sk] = pdeval(0, xmesh, psi(k, :), xmesh(end));
    psi_x_S(k) = psi_x_Sk / s_curr;
    if k > 1
        psi_t_S(k) = (psi_S(k) - psi_S(k-1))/tau;
    end
end
psi_t_S(1) = psi_t_S(2);

% TODO: Compute these using second derivative information through PDE.
% This would allow Adjoint and Forward to use nearly the same post-processing
% step.
% In a similar way to above, we would like to produce values of psi_t(x,t),
% but pdeSolver produces values \tilde{psi}(y,t)=psi(ys(t),t)
% Hence we return
% psi_t(x,t) = - x s'(t) / s^2(t) \tilde{\psi}_y(y,t) + \tilde{\psi}_t(y,t)
%            = - y s'(t)/s(t) \tilde{\psi}_y(y,t) + \tilde{\psi}_t(y,t)
% We need only the trace at x=s(t), which corresponds to y=1
psi_t_S = psi_t_S - s_der .* psi_x_S;

end