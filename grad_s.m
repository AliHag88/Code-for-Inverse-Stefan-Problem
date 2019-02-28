function [grad]=grad_s(u,psi,psi_x,tmesh,s_der,svals,u_T,mu,u_x_S,psi_x_S,psi_t_S,psi_t,psi_S,u_S,au_xx_S,s_star,psi_T)

% s_{k+1}(t)= s_k(t)-alpha_k*2*[(u_k-mu)*u_kx+phi_k+phi_kx*s'_k+phi_kt-phi_k(a*u_kx)_x]|_{x=s_k(t)}
% s_{k+1}(T)=s_k(T)-alpha_*[|u_k(s_k(T),T)-omega(s_k(T))|^2+2*(s_k(T)-s_{*})-phi_k|_{(s_k(T),T)}]

 
% grad=rand(length(tmesh),1);
grad=(2*(u_S-mu).*u_x_S+ psi_x_S.*s_der+psi_t_S-psi_S.*au_xx_S);
grad(end)=(abs(u_S(end)-u_T(end)))^2+2*(svals(end)-s_star)-psi_T(end);
% grad(1)=0;
% grad(end)=0;
grad=grad';
end


