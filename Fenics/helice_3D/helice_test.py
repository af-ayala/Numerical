from fenics import *
import matplotlib.pyplot as plt

plt.figure(figsize=(10,10))
mesh = Mesh("helice.xml")
plot(mesh, color='cyan')

plt.show()