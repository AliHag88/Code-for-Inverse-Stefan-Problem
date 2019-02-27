function [u_true_T,mu,s_true_vec,s_star,s_true,u_true]=true_solution(xmesh,tmesh)
% xmesh = linspace(0,1,100);  % Space discretization for forward problem
len_xmesh = length(xmesh); % Number of space grid points
% tmesh = linspace(0,1,100); % Time discretization for forward problem
len_tmesh = length(tmesh); % Number of time grid points

% time stepsize
tau=(tmesh(2)-tmesh(1));

% true s(t)

s_true=@(t) t.^2.*sqrt(t+0.2)+1;
% s_true=@(t) 1+t.*cos(t);

% vector values of true s(t)

s_true_vec=s_true(tmesh);

% s(T) - s(t) at final time

s_star=s_true((length(tmesh)-1)*tau);

% Define vector of s values
svals = s_true(tmesh);

% Define the vector of derivatives for s values
s_der=zeros(len_tmesh, 1);
s_der(1)=(svals(2)-svals(1))/tau;
s_der(len_tmesh(1))=(svals(len_tmesh)-svals(len_tmesh-1))/tau;
 for i=2:len_tmesh-1
     s_der(i)=(svals(i)-svals(i-1))/tau;
 end
 
% Define vector of a values
avals = 1*ones(len_tmesh, 1);

% Interpolate a values to form a(t)
af = @(t) interp1(tmesh, avals, t);

% Interpolate s' values to form s(t) derivatives
sder = @(t) interp1(tmesh, s_der, t);

% Calculate solution on rectangular domain
u_true = pdeSolver(xmesh, tmesh, af, s_true,sder);

% u(x,t) visualisation on a non-rectangular domain is based on 2d interpolation.
[X, T] = meshgrid(xmesh, tmesh);

% u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
uf = @(x,t) interp2(X, T, u_true, x/s_true(t), t, 'linear', NaN);

%TODO: Vectorize evaluation of previous function on grid.
% Just requires careful interpretation of 1/boundary_values (defined below)
% as a row-wise scaling vector.

% Define new grid on state vector's "native" domain.
boundary_values = s_true(tmesh);
x_new = linspace(0, max(boundary_values), len_xmesh);

% Initialize and evaluate u_visual, the values on the "native" domain.
u_visual = zeros(len_tmesh, len_xmesh);
for i = 1:len_tmesh
    u_visual(i, :) = uf(x_new, tmesh(i));
end

% \mu(t) Phase Transition temperature
    mu=zeros(length(tmesh),1);
    for k=1:length(tmesh)
       mu(k)=u_true(k,end);
    end
    
    
      % Values Of u(x,t) at t=T
    u_true_T=zeros(length(xmesh),1);
    for k=1:length(xmesh)
       u_true_T(k)=u_true(k,end);
    end  

end
% End Main Function

% Begin Subfunctions
function u_true = pdeSolver(xmesh, tmesh, af, sf,sder)
    m=0; % Symmetry of the problem. 0=slab (rectangular)

    % PDE to be solved is
    %    c u_t = x^{-m} (x^m f(x,t,u, DuDx))_x + s
    % The output values below give the corresponding value at each point (x, t, u, u_x)
    function [c, f, s] = pde(x,t,~,DuDx)
        c = sf(t).^2;
        f = af(t)*DuDx;
        s = sf(t)*sder(t)*x*DuDx; 
    end

    % Initial condition u(x,t_0)
    function u = ic(~)
        u = 0;
    end

    % Boundary condition at x=xl and x=xr.
    % Suffix l corresponds to left-hand side of domain,
    % Suffix r corresponds to right-hand side of domain.
    % Boundary condition takes the form
    %    p(x,t,u) + q(x,t) f(x,t,u,u_x) = 0
    % where f is defined below.
    % In particular, the inputs correspond to (x,u) at x=xl and (x,u) at x=xr
    % and the time t, and the outputs correspond to (p,q) at x=xl and (p,q) at
    % x=xr
    function [pl, ql, pr, qr] = bc(~, ~, ~, ~, t)
        pl = -sf(t);
        ql = 1;
        pr = sder(t)*sf(t);
        qr = 1;
    end

    % Calculate solution evaluated at each point of xmesh and tmesh.
    sol = pdepe(m, @pde, @ic, @bc, xmesh, tmesh);

    % Extract the first solution component as u.
    % Note: This solution is defined on tmesh * xmesh.
    % TODO: Transform to function u(x,t) within this routine.
    u_true = sol(:, :, 1);
end
% End Subfunctions