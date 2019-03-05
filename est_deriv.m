function deriv = est_deriv(vals, mesh)
  % est_deriv: Calculate first order derivative estimate for vals on mesh.
  % The output is the same shape as `vals`.
  %TODO: Instead of switching on shape here, dispatch to appropriate 
  % t- or x-derivative routine.
  if ~(isvector(vals) && isvector(mesh))
      throw(MException('est_deriv:ArgumentError', 'est_deriv requires vector input'))
  end
  if isrow(vals) && isrow(mesh)
      deriv = diff(vals) ./ diff(mesh);
      deriv(end + 1) = deriv(end);
  elseif isrow(vals) && iscolumn(mesh)
      deriv = est_deriv(vals, mesh');
  elseif iscolumn(vals) && isrow(mesh)
      deriv = est_deriv(vals, mesh');
  elseif iscolumn(vals) && iscolumn(mesh)
      deriv = (vals(2:end,1) - vals(1:end-1,1)) ./ (mesh(2:end, 1) - mesh(1:end-1,1));
      deriv(end + 1, 1) = deriv(end, 1);
  end
end
