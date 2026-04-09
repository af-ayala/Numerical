# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------
# VIGA DOBLEMENTE EMPOTRADA (FIXED-FIXED)
# Deformada dibujada hacia abajo
# --------------------------------------------

# Geometría
L = 1.0          # m
b = 0.20         # m
h = 0.20         # m

# Material
E = 25e9         # Pa

# Propiedades
I = b * h**3 / 12.0

# Coordenada
x = np.linspace(0.0, L, 400)

# --------------------------------------------
# FUNCIONES PARA VIGA FIJA-FIJA CON CARGA UNIFORME
# --------------------------------------------

def deflexion(w):
    # Signo negativo para que la deformada se grafique hacia abajo
    return -(w * x**2 * (L - x)**2) / (24.0 * E * I)

def momento(w):
    # Convención usual: momento negativo en los empotramientos
    return w * (L*x/2.0 - x**2/2.0 - L**2/12.0)

# --------------------------------------------
# CASO BASE
# --------------------------------------------
w = 60.0  # N/m

y = deflexion(w)
M = momento(w)

y_max = np.min(y)                  # más negativo = hacia abajo
y_max_abs = abs(y_max)
x_ymax = x[np.argmin(y)]

print("=== RESULTADOS ===")
print(f"Deflexión máxima = {y_max:.6e} m = {y_max*1000:.6f} mm en x = {x_ymax:.3f} m")
print(f"Magnitud de la deflexión máxima = {y_max_abs:.6e} m = {y_max_abs*1000:.6f} mm")

# --------------------------------------------
# GRAFICAS
# --------------------------------------------

# Deformada
plt.figure(figsize=(9, 4.5))
plt.plot(x, y*1000.0, linewidth=2)
plt.axhline(0.0, color='black', linewidth=0.8)
plt.scatter([x_ymax], [y_max*1000.0], s=50)
plt.text(x_ymax, y_max*1000.0, f"  Dmax = {y_max*1000:.6f} mm", va='top')
plt.title("Deformada de la viga doblemente empotrada")
plt.xlabel("x (m)")
plt.ylabel("Deflexión (mm)")
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Diagrama de momentos
plt.figure(figsize=(9, 4.5))
plt.plot(x, M, linewidth=2)
plt.axhline(0.0, color='black', linewidth=0.8)
plt.fill_between(x, M, 0.0, alpha=0.2)
plt.title("Diagrama de momentos")
plt.xlabel("x (m)")
plt.ylabel("M (N*m)")
plt.grid(True, alpha=0.3)
plt.tight_layout()

# --------------------------------------------
# BARRIDO DE CARGA
# --------------------------------------------
w_vals = np.linspace(10.0, 200.0, 20)
dmax_vals = []

for w_i in w_vals:
    y_i = -(w_i * x**2 * (L - x)**2) / (24.0 * E * I)
    dmax_vals.append(abs(np.min(y_i)) * 1000.0)   # magnitud en mm

plt.figure(figsize=(9, 4.5))
plt.plot(w_vals, dmax_vals, 'o-', linewidth=2)
plt.title("Dmax vs W_dist")
plt.xlabel("W_dist (N/m)")
plt.ylabel("Dmax (mm)")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.show()
