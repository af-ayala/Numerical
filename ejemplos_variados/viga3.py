import numpy as np
import matplotlib.pyplot as plt

# ======================
# DATOS
# ======================
L = 1.0        # longitud (m)
w = -60.0      # carga distribuida (N/m) NEGATIVA (hacia abajo)
E = 25e9       # módulo de elasticidad (Pa)
b = 0.20       # base (m)
h = 0.20       # altura (m)

# ======================
# PROPIEDADES
# ======================
I = b * h**3 / 12

# ======================
# DISCRETIZACIÓN
# ======================
x = np.linspace(0, L, 200)

# ======================
# DEFORMADA
# ======================
v = (w/(24*E*I)) * (x**4 - 2*L*x**3 + L**2*x**2)

# ======================
# DEFLEXIÓN MÁXIMA
# ======================
v_max = np.min(v)              # valor más negativo
x_max = x[np.argmin(v)]        # posición del mínimo

# ======================
# RESULTADOS
# ======================
print("========== RESULTADOS ==========")
print("Momento de inercia I =", I, "m^4")
print("Deflexión máxima =", v_max, "m")
print("Deflexión máxima =", v_max*1000, "mm")
print("Ocurre en x =", x_max, "m")
print("================================")

# ======================
# GRÁFICA
# ======================
plt.figure()

# Escalamos para visualizar mejor
plt.plot(x, v*1e6, label="Deformada")

# Línea base
plt.axhline(0)

# Punto máximo
plt.scatter(x_max, v_max*1e6, label="Deflexión máxima")

plt.xlabel("x (m)")
plt.ylabel("Deflexión (micrómetros)")
plt.title("Deformada de viga doblemente empotrada")
plt.grid()
plt.legend()

plt.show()
