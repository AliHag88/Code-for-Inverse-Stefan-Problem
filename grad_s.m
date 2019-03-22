function [grad]=grad_s(tmesh,svals,w_meas,u_T,mu_meas,u_x_S,psi_x_S,psi_t_S,psi_S,u_S,au_xx_S,s_star)
  % grad_s: Calculate gradient with respect to x=s(t)
  % Input Arguments:
  %    - tmesh: Tiem grid on which functions are evaluated
  %    - svals: Vector of boundary location values
  %    - w_meas: Function, measurement u(x,T)
  %    - u_T: Function u(x,T)
  %    - mu_meas: Function, measurement u(s(t),t)
  %    - u_x_S: Vector of values u_x(s(t),t)
  %    - psi_x_S: Vector of values psi_x(s(t),t)
  %    - psi_t_S: Vector of values psi_t(s(t),t)
  %    - psi_S: Vector of values psi(s(t),t)
  %    - u_S: Vector of values u(s(t),t)
  %    - au_xx_S: Vector of values a(t) u_{xx}(s(t),t)
  %    - s_star: Measurement s(T)
  % Output arguments:
  %    - grad: Vector of values representing J_s(t)

  t_final = tmesh(end);
  
  s_der = est_deriv(svals, tmesh);
  
  %% Test code for grad generation
  % grad=rand(length(tmesh),1);

  %% This code should mirror the gradient formula for the control x=s(t)
  % J_s(t) = [2 (u-mu) u_x - psi a u_xx + psi_x s' + psi_t]_{x=s(t)}
  % J_s(T) = |u(s(T),T)-w(s(T))|^2 + 2 (s(T)-s_*) - psi(s(T),T)
  % Note in particular that all of the vectors here should be time-like (column vectors)
  grad = (...
           2*(u_S - mu_meas(tmesh)) .* u_x_S + ...
           psi_x_S .* s_der + ...
           psi_t_S - ...
           psi_S .* au_xx_S ...
         );
  
  svals_end = svals(end);
  v1 = (u_T(svals_end)-w_meas(t_final));
  grad(end) = (...
                v1*(v1-2) + ...
                2*(svals_end-s_star)
              );
  grad(1) = 0;
end
