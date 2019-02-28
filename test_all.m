function [exitcode, passes] = test
%TEST Return true if all implemented tests pass.

  exitcode = 0;
  passes = true;

  % test_Forward
  tf1 = test_Forward(10, 20);
  tf2 = test_Forward(20, 20);
  tf3 = test_Forward(40, 40);
  passes = passes && (tf1 >= tf2) && (tf2 >= tf3);
  if ~passes
    exitcode = exitcode + 1;
    return
  end

  % test_Adjoint
  test_Adjoint_tolerance = 1e-7;
  tf1 = test_Adjoint(10, 20);
  tf2 = test_Adjoint(20, 20);
  tf3 = test_Adjoint(40, 40);
  passes = passes && tf1 < test_Adjoint_tolerance ...
           && tf2 < test_Adjoint_tolerance ...
           && tf3 < test_Adjoint_tolerance;

  if ~passes
    exitcode = exitcode + 1;
    return
  end

end

