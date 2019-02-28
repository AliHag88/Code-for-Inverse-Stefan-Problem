clc; clear all; close all;

% global t_ini t_fin Nt grad L % For Preconditioning 


% Space and Time discritization for both for forward and adjoint problem
t_ini=0;  % Initial t 
t_fin=1;  % Final t
x_ini=0;  % Initial x
x_fin=1;  % Final x
Nt=100;   % Number of t grid points
Nx=100;   % Number of x grid points
xmesh = linspace(x_ini,x_fin,Nx);  % Space discretization for forward problem
len_xmesh = length(xmesh); % Number of space grid points
tmesh = linspace(t_ini,t_fin,Nt); % Time discretization for forward problem
len_tmesh = length(tmesh); % Number of time grid points

L=1; % Preconditioning 

% time stepsize
tau=(tmesh(2)-tmesh(1));

% special stepsize

h=(xmesh(2)-xmesh(1));

% INITIAL SETUP FOR TRUE SOLUTION
initial_setup;

u_T=zeros(length(xmesh),1);
u_S=zeros(length(tmesh),1);
u_true_T=zeros(length(tmesh),1);
% TRUE SOLUTION
[u_true_T,mu,s_true_vec,s_star,s_true,u_true]=true_solution(xmesh,tmesh);

% Initial data of the true solution i.e u_true(x,t)

uInitial=u_true(1,:);

% INITIAL GUESS 

s_ini = @(t) (s_true(1)-1)*t+1;

% Initialize svals and avals

svals=s_ini(tmesh);
avals=ones(length(tmesh),1);

% alpha 

alpha=10^(-2);

% Cost Functional 

J=zeros(300,1);

k=0;

% Main Optimization Loop

while abs(svals(end)-s_star)+h*norm(u_T-u_true_T)+tau*norm(u_S-mu)<1000
   k=k+1  
   [au_xx_S,svals,s_der,u_x_S,u_S,u_T,avals,sder,u]=Forward(xmesh,tmesh,svals,avals,uInitial); 
   [psi_T,psi_t,psi_t_S,psi_S,psi_x_S,psi_x,psi]=Adjoint(xmesh,tmesh,svals,avals,u_T,u_true_T,s_der,u_S,mu);
   
   grad=grad_s(u,psi,psi_x,tmesh,s_der,svals,u_T,mu,u_x_S,psi_x_S,psi_t_S,psi_t,psi_S,u_S,au_xx_S,s_star,psi_T);
   % grad=precond()'; % Preconditioning for s(t) gradient
   svals=svals-alpha*grad;
   
   
   avals=avals-0*grad_a(u,psi,tmesh);
  
   J(k)=abs(svals(end)-s_star)+h*norm(u_T-u_true_T)+tau*norm(u_S-mu);
  
   visualization(xmesh,tmesh,svals,avals,u,u_true,len_xmesh, len_tmesh,s_true_vec,s_true,k,J)
   pause(3)
end


 