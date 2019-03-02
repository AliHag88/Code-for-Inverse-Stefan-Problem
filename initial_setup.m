% TRUE SOLUTION
[u_true_T, u_true_0, u_true_S, s_star, s_true, u_true] = true_solution(tmesh);

% alpha (size of step in anti-gradient direction)
alpha=10^(-1);

% Use true solution to set measurements
mu_meas = u_true_S;
w_meas = u_true_T;

% INITIAL GUESS
s_ini = @(t) (s_true(t_final)-1)*t+1;
a_ini = @(t) ones(size(t));

% Old information (unused)
PHI=1;

g=1;

chi=0;

gamma=1;

