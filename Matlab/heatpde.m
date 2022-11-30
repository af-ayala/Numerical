% PRIMER PASO: Definir la función que retorna los valores de la EDP
function [c,f,s] = heatdpe(x,t,u,dudx)
c = 1;
f = dudx;
s = 0;
end