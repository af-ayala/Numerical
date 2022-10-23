# ================================================
# Solución de Ecuac. Poisson 2D en Fenics
# Estamos tratando de resolver:
# -Delta u - f = 0, donde f=x*y
# ================================================

from fenics import *
import matplotlib.pyplot as plt

m = 50
def main():
    
    # Definir el dominio y malla
    mesh = UnitSquareMesh(m, m)
    
    # Definir el espacio funcional de prueba (H)
    fspace = FunctionSpace(mesh, "Lagrange",  1)
    
    # Define la función u (solución) y v (función prueba)
    u = TrialFunction(fspace)
    v = TestFunction(fspace)
    u_sol = Function(fspace)

    # Función f (rhs)
    # En fencis x=x[0], y=x[1], z=x[2]
    f = Expression('x[0]*x[1]', degree = 1)

    # Decirle a Fenics la forma variacional: a(u,v) = lv(f)
    # a(u,v) - lv(f) = 0
    a = dot(grad(u), grad(v)) * dx    
    lv = f * v * dx

    # Condición de frontera (u=0 en todo el borde)
    def funcion_bool(x, on_boundary):
        return on_boundary
    condicion_borde = DirichletBC(fspace, Constant(0.0), funcion_bool)

    # Solución con elementos finitos:
    solve(a == lv, u_sol, condicion_borde)

    # Visualización
    c = plot(u_sol, mode="color")
    plt.colorbar(c)
    # plot(mesh)

    # Adding paraview output 
    vtkfile = File('poisson/solution.pvd')
    vtkfile << u_sol
    
    #  Verificación
    print( "Verificando u(0.5,0.5) = " + str( round( u_sol(0.5,0.5) , 6) ) )
    print( "Verificando u(0.23,0.15) = " + str( round( u_sol(0.23,0.15) , 6) ) )
    
    plt.show()

if __name__ == "__main__":
    main()
