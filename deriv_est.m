function deriv = deriv_est(vals, mesh)
  % deriv_est: Calculate first order derivative estimate for vals on mesh
  deriv = diff(vals) ./ diff(mesh);
  deriv(end + 1) = deriv(end);
end
