function tf = test_optimization(len_xmesh, len_tmesh)
  % test_optimization: Run integration test for optimization code.
  % Currently, just tests that we get numerical output values for our implemented
  % test problem.
  % Input Arguments:
  %    - len_xmesh: Number of space grid points. Default: 10
  %    - len_tmesh: Number of time grid points. Default: 10
  % Set default argument values
  if ~exist('len_xmesh', 'var')
    len_xmesh = 10;
  end
  if ~exist('len_tmesh', 'var')
    len_tmesh = 10;
  end
  % Tolerance argument to pass to optimization solver. Functionally irrelevant.
  tolerance = 1e-5;
  % Number of gradient iterations to take in gradient direction (only 1 since
  % we are testing only)
  num_iterations = 1;
  % Turn off visualization in optimization routine
  do_visualization = false;

  tf = false;

  try
    [J, svals, avals] = optimization(len_xmesh, len_tmesh, tolerance, num_iterations, do_visualization);
    assert(isnumeric(J));
    assert(isnumeric(svals));
    assert(isnumeric(avals));
    tf = true;
    return
  catch ME
      switch ME.identifier
          case 'est_deriv:ArgumentError'
              formatString = 'in %s (%s)';
          otherwise
              formatString = 'internal: %s (%s)';
      end
      formatString = strcat(...
          'Exception when running optimization.m. ', ...
          formatString ...
          );
      fprintf(formatString, ME.identifier, ME.message);
  end
end % test_optimization function
