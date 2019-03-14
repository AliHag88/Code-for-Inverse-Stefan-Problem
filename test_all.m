function [exitcode, passes] = test_all
  %TEST_ALL Return exitcode==0 and passes==true if all implemented tests pass.
  % Larger (negative) exitcode implies more tests passing before failure

  %% Setup
  exitcode = 0;
  passes = true;

  %% test_Forward
  tic();
  tf1 = test_Forward(2^3, 2^3);
  time1 = toc(); tic();
  tf2 = test_Forward(2^4, 2^4);
  time2 = toc(); tic();
  tf3 = test_Forward(2^5, 2^5);
  time3 = toc();

  % Report results. Order estimate should be >1, time scaling represents increase in time
  % after 2x increase in resolution
  fprintf('Forward Problem Convergence Order Estimate: %5.3f\n', log(tf3/tf2)/log(tf2/tf1));
  fprintf('Forward Problem Time Scaling: %5.3f, %5.3f\n', time3/time2, time2/time1);

  % Fail if necessary
  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if passes
      exitcode = exitcode - 1;
  else
      return
  end

  %% test_Adjoint with zero oscillation should converge
  tic();
  tf1 = test_Adjoint(2^3, 2^3);
  time1 = toc(); tic();
  tf2 = test_Adjoint(2^4, 2^4);
  time2 = toc(); tic();
  tf3 = test_Adjoint(2^5, 2^5);
  time3 = toc();

  % Report results. Order estimate should be >1. Time scaling represents increase in time
  % after 2x increase in resolution
  fprintf('Adjoint Problem Convergence Order Estimate: %5.3f\n', log(tf3/tf2)/log(tf2/tf1));
  fprintf('Adjoint Problem Time Scaling: %5.3f, %5.3f', time3/time2, time2/time1);
  % Fail if necessary
  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if passes
      exitcode = exitcode - 1;
  else
      return
  end

  %% test_Adjoint with decreasing oscillation (should give decreasing error)
  tf1 = test_Adjoint(2^4, 2^4, 1e1);
  tf2 = test_Adjoint(2^4, 2^4, 1);
  tf3 = test_Adjoint(2^4, 2^4, 1e-1);

  % Report results. These values should be > 1
  fprintf('Adjoint Stability Estimate: Reducing Boundary Declination 10x: %5.3f, %5.3f', tf1/tf2, tf2/tf3);

  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if passes
      exitcode = exitcode - 1;
  else
      return
  end

  %% test_optimization with very small problem
  tf = test_optimization(10, 10);
  passes = passes && tf;
  if passes
      exitcode = exitcode - 1;
  else
      return
  end
  
  tf = test_true_solution(10, 10);
  passes = passes && tf;
  if passes
      exitcode = exitcode - 1;
  else
      return
  end
  
  %% Last test: run checkcode on MATLAB only
  % First, check if we're on Octave and bails out. Note that we've only gotten here
  % if all previous checks passed.
  if (exist('OCTAVE_VERSION', 'builtin') > 0)
      exitcode = 0;
      return
  end
  
  % Gather list of files to check
  checkFiles = dir('*.m');
  nfiles = length(checkFiles);
  
  % Then we run each m-file through checkcode. If the file passes, we check the cyclomatic
  % complexity of the script. More complex scripts have larger values.
  for file_i = 1:nfiles
      currFile = checkFiles(file_i).name;
      info = checkcode(currFile, '-struct', '-id');
      passes = passes && isempty(info);
      if isempty(info)
          exitcode = exitcode - 1;
          info = checkcode(currFile, '-struct', '-notok');
          if ~isempty(info)
              fprintf('Ignored %d Errors in %s. ', length(info), currFile);
          end
          cyc = checkcode(currFile, '-cyc');
          fprintf('%s\n', cyc.message);
      else
          fprintf('Checkcode failed: %s. Errors: %d.\n', currFile, length(info));
      end
  end
  if ~passes
      return
  end
  
  %% Finish up
  % If we've made it here, all tests pass. POSIX convention is exitcode==0 indicates success.
  exitcode = 0;
end

