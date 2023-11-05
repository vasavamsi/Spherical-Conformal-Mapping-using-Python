import sympy as sp
from halfedge_mesh import *
import numpy as np
from utils import *

#Obtaining the mesh in halfedge data structure for barin.off file using the Halfedge library
mesh = HalfedgeMesh("./brain.off")
vertices = mesh.vertices
edge_list = mesh.edges.keys()
faces = mesh.facets

#Creating the dictionary of neighbouring vertices
vert_dict = {}

for vertex in vertices:
    neighbouring_vertices = []
    half_edge = vertex.halfedge
    current_halfedge = half_edge
    while True:
        current_halfedge = current_halfedge.next
        if current_halfedge.vertex.index not in neighbouring_vertices:
            neighbouring_vertices.append(current_halfedge.vertex.index)
        current_halfedge = current_halfedge.opposite
        if current_halfedge is None:
            break
        else:
            if current_halfedge == half_edge:
                break
    vert_dict[vertex.index] = neighbouring_vertices

print('completed creating the vertex dictionary')

X = np.array(Gauss_map(mesh)) #Obtaining the gauss map for the input mesh to initiate the Algorithm 1.
print('completed creating the gauss map dictionary')

alpha_i = 1e-2 #intiating step size (can be changed to 1e-6 for better convergence but will take longer time to complete the convergence)
prev_Tuette_energy = 100000

#Obtaining the Tuette mapping (Algorithm 1)
"""
We can place the condition on energy difference, but as tuette embedding was becoming stable for the brain object at the value of approx 46, I utilized
this condition. With this we will save the extra loops it will take to converge.
"""
while True:
    Delta_f, Normal_comp_f, Tuette_energy = disc_tuette_energy(mesh,X,vert_dict) # The function will calculate Laplacian_f, Normal component and tuette energy
    D_f = Delta_f - Normal_comp_f
    X_next = X - alpha_i*D_f

    X_next = array_normalize(X_next)
    energy_difference = abs(Tuette_energy - prev_Tuette_energy) 
    print('Tuette energy is', Tuette_energy)
    print('Tuette energy differece', energy_difference)
    if Tuette_energy < 47:
        break
    else:
        prev_Tuette_energy = Tuette_energy
        X = X_next

#Storing the tuette mapping to .off format
file_path = './results/Spherical_conformal_mapping.off'
new_updated_vertices = X.tolist()
new_faces = []
for face in faces:
    new_faces.append([face.a,face.b,face.c])
file_printer(file_path,new_updated_vertices,new_faces,obj_output=False)

#Obtaining the Conformal mapping (Algorithm 2)
prev_Harm_energy = Tuette_energy
energy_difference = 1   
h = X_next
while energy_difference > 1e-5:
    Dh = abs_derivative(mesh,h,vert_dict) 
    h_next = h - Dh*alpha_i
    h_next = h_next - C_vector(h_next,mesh)
    h_next = array_normalize(h_next)
    harm_energy = harmonic_energy(h_next,mesh,vert_dict)
    print('Harmonic_energy', harm_energy)
    energy_difference = prev_Harm_energy - harm_energy
    print(energy_difference)
    prev_Harm_energy = harm_energy
    h = h_next

#Storing the final mapping to .off format
file_path = './results/Spherical_conformal_mapping.off'
new_updated_vertices = h.tolist()
file_printer(file_path,new_updated_vertices,new_faces,obj_output=False)