function [psi_T,psi_t,psi_t_S,psi_S,psi_x_S,psi_x,psi]=Adjoint(xmesh,tmesh,svals,avals,u_T,u_true_T,s_der,u_S,mu)
% Adjoint: Compute values of adjoint for ISP
% Input Arguments:
%    - xmesh: Space discretization
%    - tmesh: Time discretization
%    - svals: Position of boundary at time grid points
%    - avals: Diffusion coefficient a(t) at time grid points
%    - u_T: Vector of values u(x,T;v)
%    - u_true_T: Vector of measurements u(x,T)=:w(x)
%    - s_der: Vector of values s'(t)
%    - u_S: Vector of trace values u(s(t),t;v)
%    - mu: Vector of values of measurement u(s(t),t)=:mu(t)

len_xmesh = length(xmesh); % Number of space grid points
len_tmesh = length(tmesh); % Number of time grid points

h=xmesh(3)-xmesh(2);
tau=tmesh(3)-tmesh(2);

svals=fliplr(svals);

s_der=fliplr(s_der);

u_S=fliplr(u_S);

mu=fliplr(mu);

% Interpolate s values to form s(t)
s_ini = @(t) interp1(tmesh, svals, t);

% Interpolate s' values to form s'(t)
sder = @(t) interp1(tmesh, s_der, t);

% Interpolate a values to form a(t)
af = @(t) interp1(tmesh, avals, t);

uT=@(x) interp1(xmesh,u_T,x);

uTrueT=@(x) interp1(xmesh,u_true_T,x);

uS=@(t) interp1(tmesh, u_S, t);

mU=@(t) interp1(tmesh, mu, t);

%%%
%%% Calculate solution on rectangular domain
%%%
psi = pdeSolver(xmesh, tmesh, af, s_ini, sder, uT, uTrueT, uS, mU);

%%%
%%% Calculate values derived from psi
%%%
% psi_x (x,t) matrix of derivatives of psi(x,t)
psi_x=zeros(len_tmesh, len_xmesh);
psi_x(:, 1:end-1) = psi(:, 2:end) - psi(:, 1:end-1);
psi_x(:, end) = psi(:, end) - psi(:, end-1);
% Scale all values
psi_x = psi_x / h;

for i=1:len_xmesh
psi_x(:,i)=psi_x(:,i)./svals';
end

% Values of psi_x (x,t) at s(t)
psi_x_S=psi_x(:, end);


% psi_t (x,t) matrix of derivatives of psi(x,t) w.r.t t
psi_t=zeros(length(tmesh),length(xmesh));
psi_t(1:end-1, :) = psi(2:end,:) - psi(1:end-1,:);
psi_t(end, :) = psi(end, :) - psi(end-1,:);
% Scale all values
psi_t = psi_t / tau;

for i=1:len_xmesh
    psi_t(:,i)=psi_t(:,i)-psi_x(:,i).*s_der;
end

% Values of psi_t (x,t) at s(t)
psi_t_S = psi_t(:, end);

% Values of psi (x,t) at s(t)
psi_S=psi(:, end);

% Values Of psi(x,t) at t=T
psi_T=psi(end, :);

end
% End Main Function


% Begin Subfunctions
function psi = pdeSolver(xmesh, tmesh, af, sf,sder,uT,uTrueT,uS,mU)
m=0; % Symmetry of the problem. 0=slab (rectangular)

% PDE to be solved is
%    c u_t = x^{-m} (x^m f(x,t,u, DuDx))_x + s
% The output values below give the corresponding value at each point (x, t, u, u_x)
    function [c, f, s] = pde(x,t,~,DuDx)
        c = -sf(t)^2;
        f = -af(t)*DuDx;
        s =sf(t)*sder(t)*x*DuDx;
    end

% Initial condition u(x,t_0)
    function psi = ic(x)
        psi = 2*(uT(x)-uTrueT(x));
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
    function [pl, ql, pr, qr] = bc(~, ul, ~, ~, t)
        pl =0;
        ql = 1;
        pr =sder(t)*sf(t)*ul+2*sf(t)*(uS(t)-mU(t));
        qr = -1;
    end

% Calculate solution evaluated at each point of xmesh and tmesh.
sol = pdepe(m, @pde, @ic, @bc, xmesh, tmesh);

% Extract the first solution component as u.
% Note: This solution is defined on tmesh * xmesh.
% TODO: Transform to function u(x,t) within this routine.
psi = flipud(sol(:, :, 1));
end
% End Subfunctions