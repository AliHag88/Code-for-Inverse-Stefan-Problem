function [grad]=grad_a(u,psi,tmesh,xmesh)%#ok<INUSL>
  
  u_x=est_x_partial(u, xmesh(2)-xmesh(1));
  u_xx=est_x_partial(u_x, xmesh(2)-xmesh(1));
  ux_x_0=u_x(:,1); 
  psi_x_0=psi(:,1);
  ux_x_S=u(:,end);
  psi_x_S=psi(:,end);
  
  grad = -ux_x_0.*psi_x_0-ux_x_S.*psi_x_S;
end


