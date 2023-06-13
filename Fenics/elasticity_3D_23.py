# %%
from fenics import *
import matplotlib.pyplot as plt

# %%
# Variables del problema
# Viga en voladizo de sección transversal 0.2 x 0.2
L = 1
W = 0.2

mu = 1
lambda_ = 1.25

rho = 1      # unidad: kg/m3
g =  0.016   # aceleracion debida a la gravedad

# Creamos una malla para la geometria de tipo prisma, discretizando en 10x3x3
mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 10, 3, 3)
plot(mesh)

# %%
# Definimos la condicion de borde en el extremo izquierdo,
V = VectorFunctionSpace(mesh, 'P', 1)

tol = 1E-14
# Sea p = (x,y,z) un punto de R^3, 
# Recordar que en fenics: x=p[0], y=p[1], z=p[2] 

# Esta función retorna VERDADERO o FALSO
def empotramiento(p, on_boundary):
    return on_boundary and p[0]<tol

# A continuación establacemos la condición de que u vale 0 en el empotramiento:
bc = DirichletBC(V, Constant((0,0,0)), empotramiento)

# %%
# Lo siguiente no cambiarlo:

# Forma variacional de la EDP: a(u,v) = Lv(f):  <σ(u), ε(v)> = <f, v> + <T, v>
u = TrialFunction(V)
v = TestFunction(V)

d = u.geometric_dimension() # En este caso 3 dimensiones

# Definir la deformacion ε:
def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)

# Definir los esfuerzos σ = λ tr(ε) I₃ + 2 μ ε
def sigma(u):
    return lambda_* tr(epsilon(u))*Identity(d) + 2*mu*epsilon(u)

# defino mi fuerza volumétrica o de masa:
f = Constant((0,0,-rho*g))

# fuerzas de superficie:
T = Constant((0,0,0))

# %%
# Forma Variacional
a = inner(sigma(u), epsilon(v))*dx
L = dot(f,v)*dx + dot(T,v)*dx

# Calculamos la solucion
u_sol = Function(V)
solve(a==L, u_sol, bc)

# %%
# plot(u_sol, title='Desplazamiento', mode='displacement')
plot(u_sol)

# %%
# Gráfico de los esfuerzos
s = sigma(u_sol) - (1./3)*tr(sigma(u_sol))*Identity(d)
von_Mises = sqrt(3./2*inner(s, s))

V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)
plot(von_Mises, title='Esfuerzos de Von Mises')

# %%
# Cálculo de la magnitud total del desplazamiento ||u||, |u|
u_magnitude = sqrt(dot(u_sol, u_sol))
u_magnitude = project(u_magnitude, V)

plot(u_magnitude, 'Magnitud de desplazamiento')
print('min/max u:',
u_magnitude.vector().min(),
u_magnitude.vector().max())

# %%
# Escribimos los resultados para visualizarlos en Paraview
u_sol.rename("Desplazamientos", "")
von_Mises.rename("Esfuerzos de Von Mises", "")

# Exportamos a Paraview
archivo_deflexion = XDMFFile("viga_3D_2023/deflexion.xdmf")
archivo_deflexion.parameters["flush_output"] = True
archivo_deflexion.parameters["functions_share_mesh"] = True

archivo_deflexion.write(u_sol, 0.0)
archivo_deflexion.write(von_Mises, 0.0)



