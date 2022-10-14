import fenics as fe
import matplotlib.pyplot as plt

N_POINTS_P_AXIS = 12
FORCING_MAGNITUDE = 1.0

def main():
    # Mesh and Finite Element Discretization
    mesh = fe.UnitSquareMesh(N_POINTS_P_AXIS, N_POINTS_P_AXIS)
    lagrange_polynomial_space_first_order = fe.FunctionSpace(
        mesh,
        "Lagrange",  # "CG"
        1,
    )

    # Boundary Conditions
    def boundary_boolean_function(x, on_boundary):
        return on_boundary

    homogeneous_dirichlet_boundary_condition = fe.DirichletBC(
        lagrange_polynomial_space_first_order,
        fe.Constant(0.0),
        boundary_boolean_function,
    )

    # Trial and Test Functions
    u_trial = fe.TrialFunction(lagrange_polynomial_space_first_order)
    v_test = fe.TestFunction(lagrange_polynomial_space_first_order)

    # Weak Form
    forcing = fe.Constant(- FORCING_MAGNITUDE)
    weak_form_lhs = fe.dot(fe.grad(u_trial), fe.grad(v_test)) * fe.dx
    weak_form_rhs = forcing * v_test * fe.dx

    # Finite Element Assembly and Linear System solve
    u_solution = fe.Function(lagrange_polynomial_space_first_order)
    fe.solve(
        weak_form_lhs == weak_form_rhs,
        u_solution,
        homogeneous_dirichlet_boundary_condition,
    )

    # Visualize
    c = fe.plot(u_solution, mode="color")
    plt.colorbar(c)
    fe.plot(mesh)

    plt.show()

if __name__ == "__main__":
    main()