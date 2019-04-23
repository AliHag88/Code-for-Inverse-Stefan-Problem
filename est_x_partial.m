function ux = est_x_partial(vals, xmesh_spacing)
  % est_x_partial: Estimate the x-derivative of the matrix vals, assuming
  % that vals(i, :) represents a space-like line in the domain, and
  % the x-mesh spacing is uniformly `xmesh_spacing`.
  % Implemented using a first-order difference at each grid point.

  %TODO: Use built-in diff routine, selecting the appropriate dimensions
  ux = zeros(size(vals));
  ux(:, 1:end-1) = vals(:, 2:end) - vals(:, 1:end-1);
  ux(:, end) = vals(:, end) - vals(:, end - 1);
  ux = ux / xmesh_spacing;
end
