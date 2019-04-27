% Finite difference method for 1D Heat equation BIVP
% Explicit scheme
% Grid end points
x_begin = 0;
x_end = 2;
t_begin = 0;
t_end = 0.1;

% Time & space step
J = 10;
N = 60;
dt = t_end/N; % time
dx = x_end/J; % space

mu = dt/(dx^2);

x = x_begin:dx:x_end;
u = zeros(J + 1, N); % Solution grid
f = 2*x.^2; % Initial condition function u(x, 0)

for n = 1:N
    t = n*dt; % Time update
    u_bc_begin = 0; % u(0, t) = 0
    u_bc_end = 8; % u(1, t) = 8
    
    u(1, n) = u_bc_begin;
    u(J + 1, n) = u_bc_end;
    if n == 1
        for j = 2:J
            u(j, n) = mu*f(j - 1) + (-2*mu + 1)*f(j) + mu*f(j + 1);
        end
    else
        for j = 2:J
            u(j, n) = mu*u(j + 1, n - 1) + (-2*mu + 1)*u(j, n - 1) + mu*u(j - 1, n - 1);
        end 
    end
end 

tt = dt:dt:N*dt;
surf(x, tt, u');
title('Finite difference method for 1D Heat equation')
xlabel('\it x', 'interpreter', 'latex')
ylabel('\it y', 'interpreter', 'latex')
zlabel('\it z', 'interpreter', 'latex')