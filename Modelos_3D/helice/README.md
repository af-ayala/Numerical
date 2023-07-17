# Desde Autocad, Civil o Similar, pasar a STL (o similar) -> .mesh

open .mesh
create physical groups for surfaces and volumes
export as .msh, ASCII 2 

then do
dolfin-convert helice.msh helice.xml

This will create:
    - helice_facet_region.xml 
    - helice_physical_region.xml  
    - helice.xml
