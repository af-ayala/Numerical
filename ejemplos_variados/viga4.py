import numpy as np
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# ======================
# DATOS FIJOS DE LA VIGA
# ======================
L = 1.0        # longitud (m)
E = 25e9       # módulo de elasticidad (Pa)
b = 0.20       # base (m)
h = 0.20       # altura (m)

# Momento de inercia
I = b * h**3 / 12

# Coordenadas a lo largo de la viga
x = np.linspace(0, L, 300)

# ======================
# FUNCIÓN DEFORMADA
# ======================
def deformada(w):
    return (w / (24 * E * I)) * (x**4 - 2 * L * x**3 + L**2 * x**2)

# ======================
# FUNCIÓN PARA ACTUALIZAR
# ======================
def actualizar_grafica():
    try:
        w = float(entry_carga.get())
    except ValueError:
        resultado_var.set("Ingresa un número válido.")
        return

    v = deformada(w)

    # Máxima deflexión:
    # si la carga es negativa, buscamos el mínimo
    # si la carga es positiva, buscamos el máximo
    if w < 0:
        idx = np.argmin(v)
    else:
        idx = np.argmax(v)

    x_max = x[idx]
    v_max = v[idx]

    # Limpiar eje
    ax.clear()

    # Graficar deformada escalada para que se vea
    ax.plot(x, v * 1e6, linewidth=2, label="Deformada")
    ax.axhline(0, linewidth=1)
    ax.scatter([x_max], [v_max * 1e6], s=50, label="Deflexión máxima")

    ax.set_xlabel("x (m)")
    ax.set_ylabel("Deflexión (micrómetros)")
    ax.set_title("Deformada de viga doblemente empotrada")
    ax.grid(True)
    ax.legend()

    canvas.draw()

    resultado_var.set(
        f"Carga = {w:.2f} N/m | Deflexión máxima = {v_max:.6e} m = {v_max*1000:.6e} mm | x = {x_max:.3f} m"
    )

# ======================
# VENTANA PRINCIPAL
# ======================
root = tk.Tk()
root.title("Viga doblemente empotrada")

# Marco superior
frame_superior = ttk.Frame(root, padding=10)
frame_superior.pack(side=tk.TOP, fill=tk.X)

ttk.Label(frame_superior, text="Coloca la carga w (N/m):").pack(side=tk.LEFT, padx=5)

entry_carga = ttk.Entry(frame_superior, width=15)
entry_carga.pack(side=tk.LEFT, padx=5)
entry_carga.insert(0, "-60")

btn_graficar = ttk.Button(frame_superior, text="Graficar", command=actualizar_grafica)
btn_graficar.pack(side=tk.LEFT, padx=5)

# Resultado
resultado_var = tk.StringVar()
resultado_var.set("Ingresa la carga y presiona Graficar.")

label_resultado = ttk.Label(root, textvariable=resultado_var, padding=10, wraplength=900)
label_resultado.pack(side=tk.TOP, fill=tk.X)

# Figura de matplotlib embebida en tkinter
fig = Figure(figsize=(8, 5), dpi=100)
ax = fig.add_subplot(111)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Permitir Enter para actualizar
root.bind("<Return>", lambda event: actualizar_grafica())

# Graficar al inicio
actualizar_grafica()

# Ejecutar ventana
root.mainloop()
