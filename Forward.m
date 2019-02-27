function [au_xx_S,svals,s_der,u_x_S,u_S,u_T,avals,sder,u] = Forward(xmesh,tmesh,svals,avals,uInitial)
% Forward: Solve forward problem on given mesh with given data.
% Input Arguments:
%    - xmesh: Space discretization on [0,1]
%    - tmesh: Time discretization
%    - svals: Boundary locations
%    - avals: Values of diffusion coefficient a(t)
%    - uInitial: Function u(x,0)
%
% Output Arguments:
%    - u: Matrix of output values. Rows are time-levels, columns are space-levels

% TODO: Remove avals, svals, and (?) sder and s_der from output arguments

len_xmesh = length(xmesh); % Number of space grid points
len_tmesh = length(tmesh); % Number of time grid points

% spatial stepsize
h = xmesh(3)-xmesh(2);

% time stepsize
tau=(tmesh(2)-tmesh(1));

% Define the vector of derivatives for s values
s_der=zeros(len_tmesh, 1);
s_der(1)=(s_der(2)-s_der(1))/tau;
s_der(len_tmesh(1))=(s_der(len_tmesh)-s_der(len_tmesh-1))/tau;
for i=2:len_tmesh-1
    s_der(i)=(s_der(i)-s_der(i-1))/tau;
end

% Interpolate s values to form s(t)
s_ini = @(t) interp1(tmesh, svals, t);

% Interpolate a values to form a(t)
af = @(t) interp1(tmesh, avals, t);

% Interpolate s' values to form s'(t)
sder = @(t) interp1(tmesh, s_der, t);

uTrue0 = @(x) interp1(xmesh, uInitial, x);

%%%
%%% Calculate solution on rectangular domain
%%%
u = pdeSolver(xmesh, tmesh, af, s_ini, sder, uTrue0);

% Define "new" (non-rectangular) boundary curve
% s_new =@(t) interp1(tmesh, svals, t);

% u(x,t) visualisation on a non-rectangular domain is based on 2d interpolation.
% [X, T] = meshgrid(xmesh, tmesh);

% u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
% uf = @(x,t) interp2(X, T, u, x/s_new(t), t, 'linear', NaN);

avals=af(tmesh');

%%%
%%% Calculate values of solution at requested locations
%%%

%% Values Of u(x,t) at t=T
u_T=u(end, :);

%% Values of u(x,t) at s(t)
u_S=u(:, end);

%% u_x (x,t) matrix of derivatives of u(x,t)
u_x=zeros(len_tmesh, len_xmesh);

% Use forward difference for internal points
u_x(:, 1:end-1) = u(:, 2:end) - u(:, 1:end-1);

% Use backward difference at the right-hand side
u_x(:, end) = u(:, end) - u(:, end - 1);

% Scale all at once
u_x = u_x / h;

for i=1:len_xmesh
u_x(:,i)=u_x(:,i)./svals';
end


% Values of u_x (x,t) at s(t)
u_x_S=u_x(:, end);

%% Values of a(t)u_x(x,t)
au_x=zeros(len_tmesh,len_xmesh);
for i=1:len_xmesh
    au_x(:,i)=avals.*u_x(:,i);
end


%% values of (a(t)u_x(x,t))_x
au_xx=zeros(len_tmesh, len_xmesh);

% Calculate values inside domain using forward difference
for i=1:len_tmesh
    for j=1:len_xmesh-1
        au_xx(i,j)=au_x(i,j+1)-au_x(i,j);
    end
end

% Calculate values on right-hand side of domain using backward difference
au_xx(:, end) = au_x(:, end) - au_x(:, end - 1);

% Scale all differences at once
au_xx = au_xx / h;


for i=1:len_xmesh
au_xx(:,i)=au_xx(:,i)./((svals').^2);
end



% values of (a(t)u_x(x,t))_x at s(t)
au_xx_S=au_xx(:, end);

end
%%% End Main Function

%%% Begin Subfunctions
function u = pdeSolver(xmesh, tmesh, af, sf, sder, uTrue0)
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
    function u = ic(x)
        u = uTrue0(x);
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
    u = sol(:, :, 1);
    
end
%%% End Subfunctions