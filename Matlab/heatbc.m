% PASO 3: DEFINIR LAS CONDICIONES DE BORDEz
function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur - 1;
qr = 0;
end