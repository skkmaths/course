%------------------------------------------------------------------------------
% Solve linear convection equation with periodic BC
% N      = number of grid points
% cfl    = cfl number
% scheme = BD, FD, LW   | BD---> FTBS, FD----> FTFS, LW ---> Lax-Wendroff
%------------------------------------------------------------------------------
solve_periodic(101,0.9, 'bd');  % Calling the main function to execute the code
                                 % By varying the arguements here we can
                                 % solve using different schemes

function solve_periodic(N, cfl, scheme)

xmin = 0;
xmax = 1;
a    = 1;
Tf   = 2;

h  = (xmax - xmin)/(N-1);
%dt = cfl * h / abs(a);
dt = 1.1 *h;
nu = a * dt / h;

fprintf(1,'N      = %d\n', N);
fprintf(1,'cfl    = %f\n', cfl);
fprintf(1,'scheme = %s\n', scheme);
fprintf(1,'h      = %f\n', h);
fprintf(1,'dt     = %f\n', dt);
fprintf(1,'nu     = %f\n', nu);

% Make grid
x  = linspace(xmin, xmax, N);

% Initial condition
f = @(x) sin(2*pi*x);
xe = linspace(xmin, xmax, 200);

% Set initial condition
u = f(x);

t = 0;
while t < Tf

    % if t +dt > 10
    % dt = 10-t
    % nu = adt/h

   if strcmp(scheme,'bd')
      u = update_bd(nu, u);
   elseif strcmp(scheme, 'fd')
      u = update_fd(nu, u);
   elseif strcmp(scheme, 'lw')
      u = update_lw(nu, u);
   elseif strcmp(scheme, 'cs')
      u = update_ftcs(nu, u);
   else
      fprintf(1,'Unknown scheme = %s\n', scheme);
      return
   end
   t = t + dt;
   fe = f(xe - a * t);
   plot(x, u, 'b-', xe, fe, 'r-', 'LineWidth', 2)
   legend('FTBS', 'Exact')
   grid on
   pause(0.1);
end
end
%------------------------------------------------------------------------------
% Backward difference in space
%------------------------------------------------------------------------------
function u = update_bd(nu, u)

uold = u;
N = length(u);

for j=2:N
   u(j) = (1-nu)*uold(j) + nu*uold(j-1);
end

u(1) = u(N); % periodic boundary condition
end
%------------------------------------------------------------------------------
% Forward difference in space
%------------------------------------------------------------------------------
function u = update_fd(nu, u)

uold = u;
N = length(u);

for j=2:N-1
   %u(j) = (1+nu)*uold(j) - nu*uold(j+1);
   u(j) = 0.5*(1-nu)*uold(j+1) + 0.5*(1+nu)*uold(j-1);
end
u(1) = 0.5*(1-nu)*uold(j+1) + 0.5*(1+nu)*uold(N-1);
u(N) = u(1);
end
%--------------------------------------
%Forward Time Centre Space
%_______________________________________
function u = update_ftcs(nu, u)
uold = u;
N= length(u);
j=1;
u(j) = uold(j)- 0.5*nu*(uold(j+1)-uold(N-1));

for j=2:N-1
   u(j) = uold(j)- 0.5*nu*(uold(j+1)-uold(j-1));
end
u(N) = u(1);
end
%------------------------------------------------------------------------------
% Lax-Wendroff
%------------------------------------------------------------------------------
function u = update_lw(nu, u)

uold = u;
N = length(u);

j = 1;
u(j) = uold(j) - 0.5*nu*(uold(j+1) - uold(N-1)) ...
     + 0.5*nu^2*(uold(N-1) - 2*uold(j) + uold(j+1));

for j=2:N-1
   u(j) = uold(j) - 0.5*nu*(uold(j+1) - uold(j-1)) ...
        + 0.5*nu^2*(uold(j-1) - 2*uold(j) + uold(j+1));
end

u(N) = u(1);
end
