#! /bin/bash
#
./dolfin-convert cylinder_2d.msh cylinder_2d.xml
./dolfin-convert cylinder_3d.msh cylinder_3d.xml
./dolfin-convert ell_2d.msh ell_2d.xml
./dolfin-convert pinch.msh pinch.xml
./dolfin-convert step_2d.msh step_2d.xml
./dolfin-convert step_3d.msh step_3d.xml
./dolfin-convert test02.msh test02.xml
./dolfin-convert test03.msh test03.xml
#
echo "Normal end of execution."
