function [exitcode, passes] = test_all
  %TEST_ALL Return exitcode==0 and passes==true if all implemented tests pass.

  %% Setup
  exitcode = 0;
  passes = true;

  %% test_Forward
  tf1 = test_Forward(10, 20);
  tf2 = test_Forward(20, 20);
  tf3 = test_Forward(40, 40);
  
  % Report results
  fprintf('Forward Problem Convergence Order Estimate %5.3f', log(tf3/tf2)/log(tf2/tf1)); 
  
  % Fail if necessary
  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if ~passes
    exitcode = exitcode + 1;
    return
  end

  %% test_Adjoint with zero oscillation should converge
  tf1 = test_Adjoint(10, 20);
  tf2 = test_Adjoint(20, 20);
  tf3 = test_Adjoint(30, 30);
  
  % Report results
  fprintf('Adjoint Problem Convergence Order Estimate %5.3f', log(tf3/tf2)/log(tf2/tf1));
  
  % Fail if necessary
  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if ~passes
    exitcode = exitcode + 1;
    return
  end

  %% test_Adjoint with decreasing oscillation (should give decreasing error)
  tf1 = test_Adjoint(20, 20, 1e1);
  tf2 = test_Adjoint(20, 20, 1);
  tf3 = test_Adjoint(20, 20, 1e-1);
  
  % Report results
  fprintf('Adjoint Stability Estimate: Reducing Boundary Declination 10x: %5.3f, %5.3f', tf1/tf2, tf2/tf3);
  
  passes = passes && (tf1 >= tf2);

end

