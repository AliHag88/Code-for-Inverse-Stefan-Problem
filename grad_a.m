function [grad]=grad_a(u,psi,tmesh,xmesh,svals)

  grad = zeros(size(tmesh));
  for ti = 1:length(tmesh)
    xmesh_curr = xmesh*svals(ti);
    [u_trace, u_x_trace] = pdeval(0, xmesh_curr, u(ti, :), xmesh_curr);
    [psi_trace, psi_x_trace] = pdeval(0, xmesh_curr, psi(ti, :), xmesh_curr);
    grad(ti) = -2*u_x_trace(1)*psi_trace(1) - trapz(xmesh_curr, u_x_trace .* psi_trace);
  end
end

