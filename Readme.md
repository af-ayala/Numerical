# Instalar Fenics en Windows - 2023

1. Installar Virtual Box 
https://www.virtualbox.org/wiki/Downloads

2. Crear su ambiente de Linux
https://www.youtube.com/watch?v=ZrHVHyIt-PY 

3. Vaya a su terminal y tipee lo siguiente

~~~
sudo add-apt-repository ppa:fenics-packages/fenics
~~~

~~~
sudo apt update
sudo apt upgrade
sudo apt-get install fenics
sudo apt --fix-broken install   
~~~

4. Pruebe:
~~~
cd Fenics/
python3 poisson_fenics.py
~~~

https://aula.upfic.pe/course/newzoom.php?sala=99913274233

# Install FreeFem Windows

1. Descargar el ejecutable y seguir los pasos de instalación
https://www.ljll.math.upmc.fr/lehyaric/ffcs/install.htm

# Install Gmsh Windows
1. Descargar el ejecutable de la página siguiente, seleccionando la opción Windows
https://gmsh.info/

# Install Paraview Windows
1. Descargar el ejecutable de la página siguiente, seleccionando la opción Windows
https://www.paraview.org/download/


# Trabajando con mallas en Fenics y Freefem

1. Ejemplo de un dique en Gmsh:
~~~
https://www.pygimli.org/_examples_auto/1_meshing/plot_cad_tutorial.html
~~~
