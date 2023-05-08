function [u] = poisson_2D_pez(x,y)

N = 100;
u = 0; 


for n = 1:N
    % Integramos en x0  
    func_interior = @(x0) 2*( sin(pi*n*x).*sin(pi*n*x0) / (pi^2 * (n^2+1)) ) .*sin(pi*x0);
    quad = integral(func_interior, 0, 1);
    u = u + quad;
end

u = sin(pi*y).*u;
return 