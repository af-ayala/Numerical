'''
Resolver la ecuacion de Elasticidad para una viga en voladizo

               .+------------------------+
             .' |                      .'|
            +---+--------------------+'  |      ↓ Gravedad
Empotramda  |   |                    |   |
            |  ,+--------------------+---+
            |.'                      | .'
            +------------------------+'

Usaremos la ecuacion de la elasticidad lineal vista en clase asumiendo
que el material solo tendra una deformacion relativamente pequenia.

---------------------------------------------------------------------

- Ecuacion de Cauchy (Elasticidad):  
    − ∇⋅σ = f

- Relaciones Constitutivas Deformacion-Esfuerzo: 
    σ = λ tr(ε) I₃ + 2 μ ε

- Donde el tensor de deformacion es dado por:
    ε = 1/2 (∇u + (∇u)ᵀ)

Detalles de las variables:
σ  :  Tensor de esfuerzos (matriz de 3x3)
f  :  Fuerza (right-hand size, rhs, vector 3D)
λ  :  Coeficiende de Lame (real)
μ  :  Coeficiende Mu de Lame (real)
ε  :  Tensor de deformaciones (matriz de 3x3)
I₃ :  Matriz Identidad de 3x3
u  :  Vector 3D de desplazamiento

∇⋅ :  Operador divergencia (envia una matriz a un vector)
tr :  Operador traza (retorna la suma de elementos de la diagonal de una matriz)
∇  :  Operador gradiente (envia un vector a una matriz)
ᵀ  :  Operador transpuesta de una matriz

---------------------------------------------------------------------

Forma debil del problema:
    <σ(u), ∇v> = <f, v> + <T, v>

T:  Es conocido como vector de traccion que nos ayuda a describir condiciones de borde de tipo Neumann
    Para el caso de viga en voladizo, T=0.

Usualmente se escribe en funcion de las deformaciones (la + usada):
    <σ(u), ε(v)> = <f, v> + <T, v>


---------------------------------------------------------------------
Una vez calculamos los desplazamientos u:
u_sol = Function(V)
solve(a==L, u_sol, bc)

Esta funcion es un campo vectorial de desplazamientos. 
1. Tensor de esfuerzos
s = σ(u) − 1/3 tr(σ(u)) I₃

2. Ahora podemos los esfuerzos Von Mises, que nos define una medida escalar (en magnitud de los esfuerzos)
 σ_M = √(3/2 s : s)

---------------------------------------------------------------------


'''

from fenics import *
import matplotlib.pyplot as plt

# Variables del problema
L = 1; W = 0.2
mu = 1
rho = 1      # unidad: kg/m3
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = 9.81


# Creamos una malla para la geometria de tipo prisma
mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 10, 3, 3)
V = VectorFunctionSpace(mesh, 'P', 1)

# Definimos la condicion de borde en el extremo izquierdo
tol = 1E-14
def empotramiento(x, on_boundary):
    return on_boundary and x[0]<tol

bc = DirichletBC(V, Constant((0,0,0)), empotramiento)

# Definir la deformacion ε:
def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)

# Definir los esfuerzos σ = λ tr(ε) I₃ + 2 μ ε
def sigma(u):
    return lambda_* tr(epsilon(u))*Identity(d) + 2*mu*epsilon(u)

# Forma variacional de la EDP: a(u,v) = Lv(f):  <σ(u), ε(v)> = <f, v> + <T, v>
u = TrialFunction(V)
d = u.geometric_dimension() # En este caso 3 dimensiones
v = TestFunction(V)
f = Constant((0,0,-rho*g))
T = Constant((0,0,0))

a = inner(sigma(u), epsilon(v))*dx
L = dot(f,v)*dx + dot(T,v)*dx


# Calculamos la solucion
u_sol = Function(V)
solve(a==L, u_sol, bc)


# Calcular los esfuerzos s = σ(u) − 1/3 tr(σ(u)) I
s = sigma(u_sol) - 1/3 * tr(sigma(u_sol)) * Identity(d)

# Calcular los esfuerzos de Von Mises σ_M = √(3/2 s : s)
von_Mises = sqrt(3/2 * inner(s,s))

# Proyectas sobre un espacio de prueba
espacio_Lagrange = FunctionSpace(mesh, 'Lagrange', 1)
von_Mises = project(von_Mises, espacio_Lagrange)


# Escribimos los resultados para visualizarlos en Paraview
u_sol.rename("Desplazamiento", "")
von_Mises.rename("Esfuerzos von mises", "")

File("viga_empotrada/u.pvd") << u_sol
File("viga_empotrada/s.pvd") << von_Mises


