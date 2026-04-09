import numpy as np
import matplotlib.pyplot as plt

# Datos
L = 1.0               # m
w = 60.0              # N/m
E = 210e9             # Pa
b = 0.2               # m
h = 0.2               # m

# Momento de inercia
I = (b * h**3) / 12

# Discretización
x = np.linspace(0, L, 100)

# Deflexión
y = (w / (24 * E * I)) * x**2 * (L - x)**2

# Deflexión máxima
delta_max = np.max(y)

print(f"Deflexión máxima: {delta_max:.3e} m")

# Gráfica
plt.figure()
plt.plot(x, y)
plt.xlabel("Longitud (m)")
plt.ylabel("Deflexión (m)")
plt.title("Deformada de viga doblemente empotrada")
plt.grid()
plt.show()
