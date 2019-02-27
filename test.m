function passes = test()
%TEST Return true if all implemented tests pass.

passes = true;

% test_Forward
tf1 = test_Forward(10, 20);
tf2 = test_Forward(20, 20);
tf3 = test_Forward(40, 40);
passes = passes && (tf1 >= tf2) && (tf2 >= tf3);

end

