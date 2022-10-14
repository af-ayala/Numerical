# Install Fenics Windows 10

1. https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&utm_medium=webreferral&utm_campaign=dd-smartbutton&utm_location=header


2. Open Docker

3. Obtain Fenics:

~~~
curl -s https://get.fenicsproject.org | bash 
~~~

4. Activate Fenics:
~~~
fenicsproject run 
~~~ 

5. Test functionality:
~~~
python3 poisson_fenics.py  
~~~

