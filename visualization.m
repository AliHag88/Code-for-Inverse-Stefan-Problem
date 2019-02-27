function visualization(xmesh,tmesh,svals,avals,u,u_true,len_xmesh, len_tmesh,s_true_vec,s_true,k,J)

tau=tmesh(2)-tmesh(1);

s_ini = @(t) (s_true(1)-1)*t+1;

% subplot(2,3,1)
% surf(xmesh, tmesh, u)
% title(' Numerical solution u_k (x,t) on square domain.')
% xlabel('Distance x')
% ylabel('Time t')
% Below, we create an "image" plot, where u-values are translated to different colors.

subplot(2,3,1)
imagesc('XData',xmesh,'YData',tmesh,'CData',u)
title(' Numerical solution u_k (x,t) on square domain.')
xlabel('Distance x')
ylabel('Time t')
colorbar
axis tight
% Define "new" (non-rectangular) boundary curve
s_new =@(t) interp1(tmesh, s_true_vec, t);

% u(x,t) visualisation on a non-rectangular domain is based on 2d interpolation.
[X, T] = meshgrid(xmesh, tmesh);

% u:[0, max(s)] x [0, 1] \to R is defined by interpolation.
uf = @(x,t) interp2(X, T, u_true, x/s_new(t), t, 'linear', NaN);

%TODO: Vectorize evaluation of previous function on grid.
% Just requires careful interpretation of 1/boundary_values (defined below)
% as a row-wise scaling vector.

% Define new grid on state vector's "native" domain.
boundary_values = s_new(tmesh);
x_new = linspace(0, max(boundary_values), len_xmesh);

% Initialize and evaluate u_visual, the values on the "native" domain.
u_visual = zeros(len_tmesh, len_xmesh);
for i = 1:len_tmesh
    u_visual(i, :) = uf(x_new, tmesh(i));
end

% "image" plot of solution on its "native" domain.

subplot(2,3,2)
imagesc('XData',x_new, 'YData',tmesh, 'CData',u_visual)
hold on
plot(s_true(tmesh),tmesh,'*','color','red')
plot(svals,tmesh,'O','color','green')
plot(s_ini(tmesh),tmesh,'--','color','white')
title(' True solution u(x,t), initial s(t), true s(t) and s_k(t).')
xlabel('Distance x')
ylabel('Time t')
colorbar
axis tight


% Surface plot of solution on "native" domain.

% subplot(2,3,4)
% surf(x_new,tmesh,u_visual)
% title('Surface graph of u_k (x,t)')
% xlabel('Distance x')
% ylabel('Time t')



subplot(2,3,4)
scatter(k,J(k),'filled','MarkerFaceColor','black')
title('Cost Functional J')
hold on
 



end