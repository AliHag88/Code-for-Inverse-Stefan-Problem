% Alternative implementation using BVP solver
% function u = pdeSolver(xmesh, tmesh, af, sf, sder, g, uTrue0)
%     tau = tmesh(2) - tmesh(1);
%     u = zeros(length(tmesh), length(xmesh));
%
%     % Save initial data
%     k = 1;
%     t_curr = tmesh(k);
%     s_curr = sf(t_curr);
%     xmesh_curr = xmesh * s_curr;
%
%     u(k, :) = arrayfun(uTrue0, xmesh_curr);
%
%     % Create interpolation of initial data
%     uPrev = @(x) interp1(xmesh_curr, u(k, :), x, 'linear', 'extrap');
%     for k = 2:length(tmesh)
%         t_curr = tmesh(k);
%         a_curr = af(t_curr);
%         s_curr = sf(t_curr);
%         g_curr = g(t_curr);
%         xmesh_curr = xmesh * s_curr;
%
%         rhs = @(x, y) [
%             y(2);
%             (y(1) - uPrev(x))/(tau * a_curr);
%             ];
%         bcs = @(ya, yb) [
%             a_curr*ya(2) - g_curr;
%             a_curr*yb(2) - sder(t_curr)
%             ];
%         guess = @(x)[
%             uPrev(x);
%             -g_curr/a_curr
%             ];
%         solinit = bvpinit(xmesh_curr, guess);
%         sol = bvp4c(rhs, bcs, solinit);
%
%         uPrev = @(x) interp1(sol.x, sol.y(1, :), x, 'linear', 'extrap');
%         u(k, :) = arrayfun(uPrev, xmesh_curr);
%     end
% end
%%% End Subfunctions
