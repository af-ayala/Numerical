# Install Fenics Windows

1. Installar Virtual Box 
https://www.virtualbox.org/wiki/Downloads

2. Crear su ambiente de Linux
https://www.youtube.com/watch?v=ZrHVHyIt-PY 

3. Descargue Anaconda para Linux, para python >= 3.9
https://www.anaconda.com/products/distribution 
Descargan el
64-Bit (x86) Installer (737 MB)

Luego ir a:
https://docs.anaconda.com/anaconda/install/linux/

~~~
bash ~/Downloads/Anaconda3-2020.05-Linux-x86_64.sh
~~~


4. Corra lo siguiente en su terminal
~~~ 
conda create -n fenicsproject -c conda-forge fenics
source activate fenicsproject
~~~ 

5. Pruebe que funcione
~~~
cd UNI_MS_817/Fenics
python3 poisson_fenics.py  
~~~