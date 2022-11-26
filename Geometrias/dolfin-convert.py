'''
# Clase semana 11
# Creando un programa para convertir de gmsh a xml (formato de Fenics)
'''

import sys
    

def error(error):
    print("Error")


def gmsh2xml ( ifilename, ofilename ):

#*****************************************************************************80
#
## gmsh2xml() converts a Gmsh msh file to Dolfin XML format.
#
#  Discussion:
#
#    This function can only handle triangles and tetrahedrons.
#
#  Modified:
#
#    18 October 2014
#
    """Convert between .gmsh v2.0 format (http://www.geuz.org/gmsh/) and .xml, 
    parser implemented as a state machine:
    
        0 = read 'MeshFormat'
        1 = read  mesh format data
        2 = read 'EndMeshFormat'
        3 = read 'Nodes'
        4 = read  number of vertices
        5 = read  vertices
        6 = read 'EndNodes'
        7 = read 'Elements'
        8 = read  number of cells
        9 = read  cells
        10 = done
       
    """

    print(" Converting from Gmsh format (.msh, .gmsh) to DOLFIN XML format" )
    print(" Hello" )
#
#  Open files
#
    ifile = open(ifilename, "r")
    ofile = open(ofilename, "w")
#
#  Scan file for cell type
#
    cell_type = None
    dim = 0
    line = ifile.readline()
    while line:

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        # Read dimension
        if line.find("$Elements") == 0:
            line = ifile.readline()
            print("debug \n")
            print(line)
            num_cells  = int(line)
            num_cells_counted = 0
            if num_cells == 0:
                error("No cells found in gmsh file.")
            line = ifile.readline()
#
#  Now iterate through elements to find largest dimension.  
#
#  Gmsh format might include elements of lower dimensions in the element list.
#
#  We also need to count number of elements of correct dimensions.
#
#  Also determine which vertices are not used.
#
            dim_2_count = 0
            dim_3_count = 0
            vertices_2_used = []
            vertices_3_used = []
            while line.find("$EndElements") == -1:
                element = line.split()
                elem_type = int(element[1])
                num_tags = int(element[2])
                if elem_type == 2:
                    if dim < 2:
                        cell_type = "triangle"
                        dim = 2
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_2_used.extend(node_num_list)
                    dim_2_count += 1
                elif elem_type == 4:
                    if dim < 3:
                        cell_type = "tetrahedron"
                        dim = 3
                        vertices_2_used = None
                    node_num_list = [int(node) for node in element[3 + num_tags:]]
                    vertices_3_used.extend(node_num_list)
                    dim_3_count += 1
                line = ifile.readline()
        else:
            # Read next line
            line = ifile.readline()
#
#  Check that we got the cell type and set num_cells_counted
#
    if cell_type == None:
        error("Unable to find cell type.")

    if dim == 3:
        num_cells_counted = dim_3_count
        vertex_set = set ( vertices_3_used )
        vertices_3_used = None
    elif dim == 2:
        num_cells_counted = dim_2_count
        vertex_set = set ( vertices_2_used )
        vertices_2_used = None

    vertex_dict = {}

    for n,v in enumerate(vertex_set):
        vertex_dict[v] = n
#
# Step to beginning of file
#
    ifile.seek(0)
#
# Write header
#
    write_header_mesh(ofile, cell_type, dim)
#   
# Initialize node list (gmsh does not export all vertexes in order)
#
    nodelist = {}
#   
# Current state
#
    state = 0
#   
# Write data
#
    num_vertices_read = 0
    num_cells_read = 0

    while state != 10:

        # Read next line
        line = ifile.readline()
        if not line: break

        # Skip comments
        if line[0] == '#':
            continue

        # Remove newline
        if line[-1] == "\n":
            line = line[:-1]

        if state == 0:
            if line == "$MeshFormat":
                state = 1
        elif state == 1:
            (version, file_type, data_size) = line.split()
            state = 2
        elif state == 2:
            if line == "$EndMeshFormat":
                state = 3
        elif state == 3:
            if line == "$Nodes":
                state = 4
        elif state == 4:
            #num_vertices = int(line)
            num_vertices = len(vertex_dict)
            write_header_vertices(ofile, num_vertices)
            state = 5
        elif state == 5:
            (node_no, x, y, z) = line.split()
            if vertex_dict.has_key(int(node_no)):
                node_no = vertex_dict[int(node_no)]
            else:
                continue
            nodelist[int(node_no)] = num_vertices_read
            x = float(x)
            y = float(y)
            z = float(z)
            write_vertex(ofile, num_vertices_read, x, y, z)
            num_vertices_read +=1
            
            if num_vertices == num_vertices_read:
                write_footer_vertices(ofile)                
                state = 6
        elif state == 6:
            if line == "$EndNodes":
                state = 7
        elif state == 7:
            if line == "$Elements":
                state = 8
        elif state == 8:
            write_header_cells(ofile,  num_cells_counted)   
            state = 9
        elif state == 9:
            element = line.split()
            elem_type = int(element[1])
            num_tags  = int(element[2])
            if elem_type == 2 and dim == 2:
                node_num_list = [vertex_dict[int(node)] for node in element[3 + num_tags:]]
                for node in node_num_list:
                    if not node in nodelist:
                        error("Vertex %d of triangle %d not previously defined." %
                              (node, num_cells_read))
                n0 = nodelist[node_num_list[0]]
                n1 = nodelist[node_num_list[1]]
                n2 = nodelist[node_num_list[2]]
                write_cell_triangle(ofile, num_cells_read, n0, n1, n2)
                num_cells_read +=1 
            elif elem_type == 4 and dim == 3:
                node_num_list = [vertex_dict[int(node)] for node in element[3 + num_tags:9]]
                for node in node_num_list:
                    if not node in nodelist:
                        error("Vertex %d of tetrahedron %d not previously defined." % 
                              (node, num_cells_read))
                n0 = nodelist[node_num_list[0]]
                n1 = nodelist[node_num_list[1]]
                n2 = nodelist[node_num_list[2]]
                n3 = nodelist[node_num_list[3]]
                write_cell_tetrahedron(ofile, num_cells_read, n0, n1, n2, n3)
                num_cells_read +=1 

            if num_cells_counted == num_cells_read:
              write_footer_cells(ofile)                
              state = 10
        elif state == 10:
            break
#
# Check that we got all data
#
    if state == 10:
        print(" Conversion done" )
    else:
        error("Missing data, unable to convert \n\ Did you use version 2.0 of the gmsh file format?")
#  
# Write footer
#
    write_footer_mesh(ofile)  
#
# Close files
#
    ifile.close()
    ofile.close()


def write_header_mesh(ofile, cell_type, dim):

#*****************************************************************************80
#
## write_header_mesh() writes the mesh header.

    ofile.write("""\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>

<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">
  <mesh celltype="%s" dim="%d">
""" % (cell_type, dim))

def write_header_graph(ofile, graph_type):

#*****************************************************************************80
#
## write_header_graph() writes the graph header.
#
    ofile.write("""\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>

<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">
  <graph type="%s">
""" % (graph_type))

def write_footer_mesh(ofile):

#*****************************************************************************80
#
## write_footer_mesh() writes the mesh footer.
#
    ofile.write("""\
  </mesh>
</dolfin>
""")

def write_footer_graph(ofile):

#*****************************************************************************80
#
## write_footer_graph() writes the graph footer.
#
    ofile.write("""\
  </graph>
</dolfin>
""")

def write_header_vertices(ofile, num_vertices):

#*****************************************************************************80
#
## write_header_vertices() writes the header for the vertex information.
#
    "Write vertices header"
    print(" Expecting %d vertices" % num_vertices )
    ofile.write("    <vertices size=\"%d\">\n" % num_vertices)

def write_footer_vertices(ofile):

#*****************************************************************************80
#
## write_footer_vertices() writes the footer for the vertex information.
#
    "Write vertices footer"
    ofile.write("    </vertices>\n")
    print(" Found all vertices" )

def write_header_edges(ofile, num_edges):

#*****************************************************************************80
#
## write_header_edges() writes the header for the edge information.
#
    "Write edges header"
    print(" Expecting %d edges" % num_edges )
    ofile.write("    <edges size=\"%d\">\n" % num_edges)

def write_footer_edges(ofile):

#*****************************************************************************80
#
## write_footer_edges() writes the footer for the edge information.
#
    "Write edges footer"
    ofile.write("    </edges>\n")
    print(" Found all edges" )

def write_vertex(ofile, vertex, x, y, z):

#*****************************************************************************80
#
## write_vertex() writes the vertex information.
#
    "Write vertex"
    ofile.write("      <vertex index=\"%d\" x=\"%g\" y=\"%g\" z=\"%g\"/>\n" % \
        (vertex, x, y, z))

def write_graph_vertex(ofile, vertex, num_edges, weight = 1):

#*****************************************************************************80
#
## write_graph_vertex() ???
#
    "Write graph vertex"
    ofile.write("      <vertex index=\"%d\" num_edges=\"%d\" weight=\"%d\"/>\n" % \
        (vertex, num_edges, weight))

def write_graph_edge(ofile, v1, v2, weight = 1):

#*****************************************************************************80
#
## write_graph_edge() ???
#
	 "Write graph edge"
	 ofile.write("      <edge v1=\"%d\" v2=\"%d\" weight=\"%d\"/>\n" % \
        (v1, v2, weight))

def write_header_cells(ofile, num_cells):

#*****************************************************************************80
#
## write_header_cells() ???
#
    "Write cells header"
    ofile.write("    <cells size=\"%d\">\n" % num_cells)
    print(" Expecting %d cells" % num_cells )

def write_footer_cells(ofile):

#*****************************************************************************80
#
## write_footer_cells() ???
#
    "Write cells footer"
    ofile.write("    </cells>\n")
    print(" Found all cells" )

def write_cell_interval ( ofile, cell, n0, n1 ):
#*****************************************************************************80
#
## write_cell_interval() ???
#
    "Write cell (interval)"
    ofile.write("      <interval index=\"%d\" v0=\"%d\" v1=\"%d\"/>\n" % \
        (cell, n0, n1))

def write_cell_triangle(ofile, cell, n0, n1, n2):
#*****************************************************************************80
#
## write_cell_triangle() ???
#
    "Write cell (triangle)"
    ofile.write("      <triangle index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\"/>\n" % \
        (cell, n0, n1, n2))

def write_cell_tetrahedron(ofile, cell, n0, n1, n2, n3):

#*****************************************************************************80
#
## write_cell_tetrahedron() ???
#
    "Write cell (tetrahedron)"
    ofile.write("      <tetrahedron index=\"%d\" v0=\"%d\" v1=\"%d\" v2=\"%d\" v3=\"%d\"/>\n" % \
        (cell, n0, n1, n2, n3))

if __name__ == "__main__":
    gmsh2xml ("cube.msh", "cube.xml")

