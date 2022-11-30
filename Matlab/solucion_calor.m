% Este archivo NO necesita ser editado.
% Sólo editar los archivos heatpde.m heatic.m heatbc.m

m = 0 ; % cambie m como sea necesario.
L = 1; % la longitud del elemento de análisis
T = 5; % tiempo de simulación

% definir la malla, si desea agregue más puntos al emallado
x = linspace(0,L,100);
t = linspace(0,T,500);


sol = pdepe(m,@heatpde,@heatic,@heatbc,x,t);

colormap hot
imagesc(x,t,sol)
colorbar
xlabel('Distancia x','interpreter','latex')
ylabel('Tiempo t','interpreter','latex')
title('Ecuación de calor para 0<x<1 y $0<t<5')