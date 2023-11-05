import math
import os
from halfedge_mesh import *
import numpy as np

def midpoint(a,b,mesh): 
    """
    a,b are the vertices of the edge for which mid point is to be calculated
    mesh is the mesh object obtained using the python package provided

    returns the midpoint in array format with len 3 (x,y,z co-ordinates)
    """
    vertices = mesh.vertices
    edges = mesh.edges

    # calculating mid point for a,b edge
    if (b,a) in edges.keys():
        p1 = vertices[a].get_vertex()
        p2 = vertices[b].get_vertex()
        p3 = edges[(a,b)].prev.vertex.get_vertex()
        p4 = edges[(b,a)].prev.vertex.get_vertex()
        mid_ab = []
        mid_ab.append(0.375*(p1[0]+p2[0]) + 0.125*(p3[0]+p4[0]))
        mid_ab.append(0.375*(p1[1]+p2[1]) + 0.125*(p3[1]+p4[1]))
        mid_ab.append(0.375*(p1[2]+p2[2]) + 0.125*(p3[2]+p4[2]))
    else:
        p1 = vertices[a].get_vertex()
        p2 = vertices[b].get_vertex()
        mid_ab = []
        mid_ab.append(0.5*(p1[0]+p2[0]))
        mid_ab.append(0.5*(p1[1]+p2[1]))
        mid_ab.append(0.5*(p1[2]+p2[2]))

    return mid_ab

def updated_vertex(vertex,mesh):
    """
    vertex is the vertex of mesh for which the updated version is to be calculated
    mesh is the mesh object obtained using the python package provided

    return updated vertex of the provided original vertex
    """
    vertices = mesh.vertices
    edge_list = mesh.edges.keys()
    neighbouring_vertices_indices = []
    for i in range(len(vertices)):
        if (vertex, i) in edge_list or (i,vertex) in edge_list:
            if i not in neighbouring_vertices_indices:
                neighbouring_vertices_indices.append(i)
    n = len(neighbouring_vertices_indices)
    neighbouring_vertices = [vertices[i].get_vertex() for i in neighbouring_vertices_indices]

    ## Below commented code is written to reduce the computation due to the above loop, although it is giving the correct number of vertices and faces
    ## proper output is not obtained in the meshlab viewer.
    # neighbouring_vertices = []
    # half_edge = vertices[vertex].halfedge
    # current_halfedge = half_edge
    # neighbouring_vertices.append(current_halfedge.vertex.get_vertex())
    # while True:
    #     current_halfedge = current_halfedge.next
    #     current_halfedge = current_halfedge.opposite
    #     if current_halfedge is None:
    #         break
    #     else:
    #         if current_halfedge == half_edge:
    #             break
    #         else:
    #             neighbouring_vertices.append(current_halfedge.vertex.get_vertex())
    # # print(neighbouring_vertices)
    # n = len(neighbouring_vertices)

    #Calculating the updated vertex co-ordinates of the old vertex
    theta = math.radians(360 / n)
    cosine_part = 0.25*(math.cos(theta))
    sq_part = ( 0.375 + cosine_part)**2
    alpha = (0.625 - sq_part) / n
    V_i = vertices[vertex].get_vertex()
    # Part-1
    v_part_1 = []
    v_part_1.append((1-(n*alpha))*(V_i[0]))
    v_part_1.append((1-(n*alpha))*(V_i[1]))
    v_part_1.append((1-(n*alpha))*(V_i[2]))

    # # Part-2
    v_part_2 = []
    for i in range(3):
        tmp = 0
        for v in neighbouring_vertices:
            tmp = tmp + v[i]
        v_part_2.append(tmp)
    
    v_part_2 = [num * alpha for num in v_part_2]

    # Final updated vertex
    v_updated = []
    v_updated.append(v_part_1[0] + v_part_2[0])
    v_updated.append(v_part_1[1] + v_part_2[1])
    v_updated.append(v_part_1[2] + v_part_2[2])

    return v_updated

def file_printer(file_path,new_updated_vertices,new_faces,obj_output=False):
    """
    creates the output file in obj or off format 
    file_path = output file path in .obj or .off format
    new_updated_vertices = list of new(midpoints) & updated vertices
    new_faces = list of new face (list of list of 3 vertex indices corresponding to the new faces)
    obj_output = Boolean to decide the format of output (.obj by default)

    """

    if obj_output:

        if os.path.exists(file_path):
            print("The file '{}' already exists.".format(file_path))
        else:
            try:
                with open(file_path, 'w') as file:
                    print("File '{}' created successfully.".format(file_path))
            except IOError as e:
                print("Error: Unable to create the file '{}'.".format(str(e)))

        with open(file_path, 'w') as file:
            file.write('')

        for v in new_updated_vertices:
            with open(file_path, 'a') as file:
                file.write('v {} {} {}\n'.format(v[0]+1, v[1]+1, v[2]+1))

        for f in new_faces:
            with open(file_path, 'a') as file:
                file.write('f {}/{} {}/{} {}/{}\n'.format(f[0]+1,f[0]+1,f[1]+1,f[1]+1,f[2]+1,f[2]+1))
    else:

        if os.path.exists(file_path):
            print("The file '{}' already exists.".format(file_path))
        else:
            try:
                with open(file_path, 'w') as file:
                    print("File '{}' created successfully.".format(file_path))
            except IOError as e:
                print("Error: Unable to create the file '{}'.".format(str(e)))

        with open(file_path, 'w') as file:
            file.write('')

        with open(file_path, 'a') as file:
            file.write('OFF\n')
            file.write('{} {} 0\n'.format(len(new_updated_vertices), len(new_faces)))

        for v in new_updated_vertices:
            with open(file_path, 'a') as file:
                file.write('{} {} {}\n'.format(v[0], v[1], v[2]))
        
        for f in new_faces:
            with open(file_path, 'a') as file:
                file.write('3 {} {} {}\n'.format(f[0],f[1],f[2]))
    
def Gauss_map(mesh):
    """
    Takes in the mesh and calculates the Gauss map
    """   
    vertices = mesh.vertices
    facets = mesh.facets 
    gauss_map = []
    for vertex in vertices:
        normal = [0,0,0]
        faces_w_vertex = []
        for facet in facets:
            facet_vertices = [facet.a,facet.b,facet.c]
            if vertex.index in facet_vertices:
                faces_w_vertex.append(facet)
                facet_normal = facet.get_normal()
                # print(facet_normal)
                normal[0] = normal[0] + facet_normal[0]
                normal[1] = normal[1] + facet_normal[1]
                normal[2] = normal[2] + facet_normal[2]
        normal[0] = float(normal[0])/float(len(faces_w_vertex))
        normal[1] = float(normal[1])/float(len(faces_w_vertex))
        normal[2] = float(normal[2])/float(len(faces_w_vertex))
        normal = normalize(normal)
        gauss_map.append(normal)
    return gauss_map
    
def disc_tuette_energy(mesh, X, vert_dict):
    """
    takes in mesh, X and vert_dict to caculate the D_f and Tuette Energy
    """
    vertices = mesh.vertices
    
    # Calculating discrete laplacian and Tuette/Harmonic Energy
    Normal_comp_f = []
    Delta_f = []
    Tuette_energy = 0

    for vertex in vertices:
        delta_f = np.array([0,0,0]) #discrete laplacian
        neighbouring_vertices = vert_dict[vertex.index]
        for neighbour in neighbouring_vertices:
            k_uv = 1.0
            #Calculating Tuette energy
            vec_uv = np.array(X[vertex.index]) - np.array(X[neighbour])
            Tuette_energy = Tuette_energy + k_uv*(sum(component ** 2 for component in vec_uv))
            
            #Calculating delta_f
            delta_f = delta_f + (k_uv * (np.array(X[vertex.index]) - np.array(X[neighbour])))
        normal_comp_f = dot(delta_f,X[vertex.index])*X[vertex.index]
        Delta_f.append(delta_f)
        Normal_comp_f.append(normal_comp_f)
        # print(vertex.index)
    return(np.array(Delta_f), np.array(Normal_comp_f), Tuette_energy*0.5)

def abs_derivative(mesh, X_next,vert_dict):
    """
    Takes the mesh, vertices dictonary(keys : vertex index, values : vertex neighbours) and current mapping in the iteration to calculate the absolute
    derivative
    """
    vertices = mesh.vertices
    Delta_f = []
    for vertex in vertices:
        delta_f = np.array([0,0,0])
        f_u = X_next[vertex.index]
        for neighbour in vert_dict[vertex.index]:
            f_v = X_next[neighbour]
            #calculating kuv
            thishalfedge = vertex.halfedge
            v1 = vertex.index
            v2 = thishalfedge.prev.vertex.index
            v3 = thishalfedge.next.vertex.index

            vA = X_next[v1,:] - X_next[v3,:]
            vB = X_next[v2,:] - X_next[v3,:]
            dot_product = np.dot(vA, vB)
            magnitude1 = np.linalg.norm(vA)
            magnitude2 = np.linalg.norm(vB)
            beta = dot_product / (magnitude1 * magnitude2)
            beta = np.arccos(beta)

            opphalfedge = thishalfedge.opposite
            v4 = opphalfedge.next.vertex.index
            
            vA = X_next[v4,:] - X_next[v1,:]
            vB = X_next[v4,:] - X_next[v2,:]
            dot_product = np.dot(vA, vB)
            magnitude1 = np.linalg.norm(vA)
            magnitude2 = np.linalg.norm(vB)
            alpha = dot_product / (magnitude1 * magnitude2)
            alpha = np.arccos(alpha)

            kuv = 0.5 * ((1 / math.tan(alpha)) + (1 / math.tan(beta)))

            delta_f = delta_f + kuv*(f_v - f_u)
        Delta_f.append(delta_f)
    Delta_f = np.array(Delta_f)
    normal_comp = Delta_f*X_next
    normal_comp = np.sum(normal_comp, axis=1)
    
    normal_comp = normal_comp.reshape(-1, 1)
    # normal_comp = normal_comp * np.ones((1, len(vertices)))
    normal_comp = normal_comp*X_next
    absolute_derivative = Delta_f - normal_comp
    return absolute_derivative

def harmonic_energy(X_next,mesh,vert_dict):
    """
    Takes the mesh, vertices dictonary(keys : vertex index, values : vertex neighbours) and current mapping in the iteration to calculate the harmonic energy
    """
    vertices = mesh.vertices
    Harm_energy = 0
    for vertex in vertices:
        energy_u = 0
        f_u = X_next[vertex.index]
        for neighbour in vert_dict[vertex.index]:
            f_v = X_next[neighbour]
            #calculating kuv
            thishalfedge = vertex.halfedge
            v1 = vertex.index
            v2 = thishalfedge.prev.vertex.index
            v3 = thishalfedge.next.vertex.index

            vA = X_next[v1,:] - X_next[v3,:]
            vB = X_next[v2,:] - X_next[v3,:]
            dot_product = np.dot(vA, vB)
            magnitude1 = np.linalg.norm(vA)
            magnitude2 = np.linalg.norm(vB)
            beta = dot_product / (magnitude1 * magnitude2)
            beta = np.arccos(beta)

            opphalfedge = thishalfedge.opposite
            v4 = opphalfedge.next.vertex.index
            
            vA = X_next[v4,:] - X_next[v1,:]
            vB = X_next[v4,:] - X_next[v2,:]
            dot_product = np.dot(vA, vB)
            magnitude1 = np.linalg.norm(vA)
            magnitude2 = np.linalg.norm(vB)
            alpha = dot_product / (magnitude1 * magnitude2)
            alpha = np.arccos(alpha)

            kuv = 0.5 * ((1 / math.tan(alpha)) + (1 / math.tan(beta)))

            v1v2 = X_next[v1,:] - X_next[v2,:]
            magnitude = (np.linalg.norm(v1v2))**2
            energy_uv = kuv*magnitude
            energy_u = energy_u + energy_uv
        Harm_energy = Harm_energy + energy_u
    return (Harm_energy/2.0)

def array_normalize(X):
    """
    Takes the array of co-ordinates and returns the normalized version
    """
    X = X.tolist()
    norm_X = []
    for x in X:
        norm_X.append(normalize(x))
    return(np.array(norm_X))

def calculate_area_of_triangle(p1, p2, p3):
    """
    Takes in 3 co-ordinates of the triangle to calculate the area
    """
    # Convert the coordinates to numpy arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # Calculate two vectors in the plane of the triangle
    v1 = p2 - p1
    v2 = p3 - p1

    # Calculate the cross product of the two vectors
    cross_product = np.cross(v1, v2)

    # Calculate the magnitude of the cross product
    area = 0.5 * np.linalg.norm(cross_product)

    return area

def C_vector(X_next, mesh):
    """
    Takes in the current updated step and vertex neighbours dictonary to return the center vector for mobius transformation
    """
    vertices = mesh.vertices
    facets = mesh.facets 
    C_vec = np.array([0,0,0])
    for vertex in vertices:
        vertex_area = 0
        for facet in facets:
            facet_vertices = [facet.a,facet.b,facet.c]
            if vertex.index in facet_vertices:
                area = calculate_area_of_triangle(X_next[facet.a,:],X_next[facet.b,:],X_next[facet.c,:])
                vertex_area = vertex_area + area
        
        C_vec = C_vec + ((0.33)*vertex_area)*X_next[vertex.index,:]
    C_vector = np.full((len(vertices),3), C_vec)
    return C_vector





        
