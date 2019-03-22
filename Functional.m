function [au_xx_S, u_x_S, u_S, u_T, utilde, J] = Functional(xmesh, tmesh, svals, avals, g, u_true_0, s_star, w_meas, mu_meas)
  % Functional: Calculate functional and solution of forward problem
  % Input Arguments:
  %    - xmesh
  %    - tmesh
  %    - svals
  %    - avals
  %    - g
  %    - u_true_0
  [au_xx_S, u_x_S, u_S, u_T, utilde] = ...
      Forward(xmesh, tmesh, svals, avals, g, u_true_0);
  
  % Calculate functional
  J = abs(svals(end) - s_star)^2 + ...
        trapz(xmesh, (u_T(xmesh) - w_meas(xmesh)).^2) + ...
        trapz(tmesh, (u_S - mu_meas(tmesh)).^2);
end