% TODO: Integrate an appropriate subset of this in
% initial_setup or a similar script.
% Set default arguments
if ~exist('len_xmesh', 'var')
  len_xmesh = 100;
end
if ~exist('len_tmesh', 'var')
  len_tmesh = 40;
end
if ~exist('tolerance', 'var')
  tolerance = 1e-9;
end

% Flags below change the optimization routine
if ~exist('num_iterations', 'var')
  num_iterations = 10;
end
if ~exist('num_sub_iterations', 'var')
  num_sub_iterations = 6;
end
if ~exist('do_visualization', 'var')
  do_visualization = false;
end

% See initial_setup.m for the following parameters
if ~exist('use_synthetic_data', 'var')
  use_synthetic_data = true;
end
if ~exist('initial_data_parameter_s', 'var')
  initial_data_parameter_s = 0.6;
end
if ~exist('initial_data_parameter_a', 'var')
  initial_data_parameter_a = 0.6;
end

% Sobolev preconditioning parameters
if ~exist('sobolev_preconditioning_s', 'var')
  sobolev_preconditioning_s = 0.05;
end
if ~exist('sobolev_preconditioning_a', 'var')
  sobolev_preconditioning_a = 0.08;
end

% Choose which controls to reconstruct during optimization process
if ~exist('reconstruct_a', 'var')
  reconstruct_a = 1;
end
if ~exist('reconstruct_s', 'var')
  reconstruct_s = 1;
end

% Set final moment (as in optimization.m code)
t_final = 1;


% Thresholds for minimum values of s(t) and a(t)
svals_minimum_threshold = 1e-4;
avals_minimum_threshold = 1e-4;

% Preconditioning "mode".
% Select == 1 to impose homogeneous Neumann condition (the "original" flavor)
% Select == 2 to give a Neumann condition on the right-hand side and a
%             homogeneous Dirichlet condition on the left-hand side
s_precond_mode = 2;
a_precond_mode = 2;
