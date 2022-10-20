function [u] = poisson_2D(x,y)

N = 50;
u = 0; 

for k = 1:N
    for l = 1:N
        u = u + 4*(-1)^(k+l)*sin(k*pi*x).*sin(l*pi*y) / (k*l*(k^2 + l^2)*pi^4);
    end
end
