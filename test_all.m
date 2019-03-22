function [exitcode, passes] = test_all
  %TEST_ALL Return exitcode==0 and passes==true if all implemented tests pass.
  % Larger (negative) exitcode implies more tests passing before failure

  %% Setup
  exitcode = 0;
  passes = true;

  %% test_Forward
  tic();
  [u_err1, ~, ~] = test_Forward(2^3, 2^3);
  time1 = toc(); tic();
  [u_err2, ~, ~] = test_Forward(2^4, 2^4);
  time2 = toc(); tic();
  [u_err3, ~, ~] = test_Forward(2^5, 2^5);
  time3 = toc();

  % Report results. Order estimate should be >1, time scaling represents increase in time
  % after 2x increase in resolution
  order_est = convergence_order_est(u_err1, u_err2, u_err3);
  fprintf('Forward Problem Convergence Order Estimate: p=%5.3f\n', order_est);
  fprintf('Forward Problem Time Scaling: %5.3f, %5.3f\n', time3/time2, time2/time1);

  % Fail if necessary
  passes = passes && (u_err1 >= u_err2) && (u_err2 >= u_err3);
  if passes
      exitcode = exitcode - 1;
  else
      return
  end

  %% test_Adjoint with zero oscillation should converge
  tic();
  psi_err1 = test_Adjoint(2^3, 2^3);
  time1 = toc(); tic();
  psi_err2 = test_Adjoint(2^4, 2^4);
  time2 = toc(); tic();
  psi_err3 = test_Adjoint(2^5, 2^5);
  time3 = toc();

  % Report results. Order estimate should be >1. Time scaling represents increase in time
  % after 2x increase in resolution
  order_est = convergence_order_est(psi_err1, psi_err2, psi_err3);
  fprintf('Adjoint Problem Convergence Order Estimate: p=%5.3f\n', order_est);
  fprintf('Adjoint Problem Time Scaling: %5.3f, %5.3f', time3/time2, time2/time1);
  % Fail if necessary
  passes = passes && (psi_err1 >= psi_err2) && (psi_err2 >= psi_err3);
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
  
  %% Run some sanity checks for true solution
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

% Estimate order of convergence from e_n, e_{n+1}, and e_{n+2}.
% Note that order of convergence p and the error estimate mu are defined as
% lim_{n\to\infty} |e_{n+1}| / |e_n|^p = mu
% Hence, log(e_{n+1}) - p*log(e_n) \approx log(mu) \approx log(e_{n+2}) - p*log(e_{n+1})
% which implies log(e_{n+1}) - log(e_{n+2}) = p*(log(e_n) - log(e_{n+1}))
% and hence p = log(e_{n+1}/e_{n+2}) / log(e_n/e_{n+1})
function order_est = convergence_order_est(en, enp1, enp2)
  if abs(en) < 1e-10 || abs(enp1) < 1e-10 || abs(enp2) < 1e-10
      order_est = Inf;
      return
  end
  order_est = log(enp1/enp2)/log(en/enp1);
end
