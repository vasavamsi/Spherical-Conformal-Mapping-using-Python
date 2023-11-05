# Spherical-Conformal-Mapping-using-Python

This work is the python implementation of math work explained in the paper "Genus Zero Surface Conformal Mapping" by Yalin Wang. I completed this work as a part of Course porject in the course "CSE 570 Advanced Computer Graphics I". The project is completed using the help of https://github.com/carlosrojas/halfedge_mesh for the halfedge data structure.

The 'main.py' is the main python file to execute the subdivision algorithm. The 'results' folder stores the outputs in either ‘.obj’ or ‘.off’ format as per the variables set in the 'main.py'. The comments.

## Algorithm

1)	Converted the input mesh object ('brain.off' here) into a gauss map.
2)	We then iterate to optimize the tuette energy (algorithm 1) to get the tuette embedding.
3)	Using the tuette embedding, we iterate over to optimize the harmonic energy with the energy difference threshold as 1e-5 (can be adjusted). The converged harmonic energy for 'brain.off' comes out to be approximately 25.
4)	Create the .obj or .off file for the final conformal mapping.

##Instructions to run the code

The package runs without any errors in python = 2 environment. Change the following variables in the main.py:
1)	mesh – Change the path provided in ‘HalfedgeMesh’ function to read the input file in ‘.off’ format. 
2)	file_path – Output path to save the output file in ‘.obj’ or ‘.off’ format.
The helper functions are defined in the ‘utils.py’ file. The comments in the functions explain the input and outputs.
The obtained result is placed in the results folder.

![spherical conformal mapping](https://github.com/vasavamsi/Spherical-Conformal-Mapping-using-Python/assets/58003228/ff89d084-0b04-476c-bb9c-61943a20a0b1)

**Note:** As we have no way to include the texture information, visualization would be hard in the case of python implemantation.
