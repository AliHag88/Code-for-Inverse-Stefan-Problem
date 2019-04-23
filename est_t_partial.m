function ut = est_t_partial(vals, tmesh_spacing)
  % est_t_partial: Estimate the t-derivative of the matrix vals, assuming
  % that vals(:, i) represents a time-like line in the domain, and
  % the t-mesh spacing is uniformly `tmesh_spacing`.
  % Implemented using a first-order difference at each grid point.

  %TODO: Use built-in diff routine, selecting the appropriate dimensions
  ut = zeros(size(vals));
  ut(1:end-1, :) = vals(2:end, :) - vals(1:end-1, :);
  ut(end, :) = vals(end, :) - vals(end - 1, :);
  ut = ut / tmesh_spacing;
end
