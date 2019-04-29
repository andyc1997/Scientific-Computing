% Finite difference method for 1D Fisher's equation BIVP
% Explicit scheme
% Grid end points
x_begin = 0;
x_end = 1;
t_begin = 0;
t_end = 1;

% Time & space step

dt = 0.01; % time
dx = 0.2; % space
J = x_end/dx;
N = t_end/dt;

mu = dt/(dx^2); % In this case, stability condition is relaxed a little, however, mu should not be greater than 0.5 by a lot.

x = x_begin:dx:x_end;
u = zeros(J + 1, N); % Solution grid
f = sin(pi * x).^2; % Initial condition function u(x, 0)

for n = 1:N
    t = n*dt; % Time update
    u_bc_begin = 0; % u(0, t) = 0
    u_bc_end = 0; % u(1, t) = 0

    u(1, n) = u_bc_begin;
    u(J + 1, n) = u_bc_end;
    if n == 1
        for j = 2:J
            u(j, n) = mu*f(j - 1) + (-2*mu + 1)*f(j) + mu*f(j + 1) + dt*(f(j) - f(j)^2);
        end
    else
        for j = 2:J
            u(j, n) = mu*u(j + 1, n - 1) + (-2*mu + 1)*u(j, n - 1) + mu*u(j - 1, n - 1) + dt*(u(j, n - 1) - u(j, n - 1)^2);
        end 
    end
end 

tt = dt:dt:N*dt;
surf(x, tt, u'); % Set 'EdgeColor' to 'none' when t_F is large
zlim auto
title('Finite difference method for 1D Fisher''s equation')
xlabel('\it x', 'interpreter', 'latex')
ylabel('\it t', 'interpreter', 'latex')
zlabel('\it u(x, t)', 'interpreter', 'latex')
