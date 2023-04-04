# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 09:54:03 2023

@author: alvar
"""

from Bio.PDB import *
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


#Reading the mmcif file and getting the residues
parser = MMCIFParser()
structure = parser.get_structure("myProtein", "1dg3.cif")
residues = structure.get_residues()


#creates a list that will store the ID of each residue, and all of its atoms coordenates
all_atoms = [] 

#creates a df that will store the c alpha coordenates
calpha_atoms = pd.DataFrame(columns=['X','Y','Z']) 

for residue in residues:
    
    if residue.resname != 'HOH': #excludes the water molecules 
        
        #Stores the residue id and prints it    
        coords = [residue._id[1]]
        print(f'Residue id: {residue._id[1]} \nAtoms type and coordinates:')
        
        #Fills the all_atoms list and prints the atoms types and coordinates
        for atom in residue:
            coords.append(atom.get_coord())
            print(atom.fullname,*atom.get_coord())
        print()
        all_atoms.append(coords)
        
        #Fills the df with c alpha coordinates
        if residue.has_id("CA"):
            calpha_atoms.loc[len(calpha_atoms.index)] = residue["CA"].get_coord()
            

#Plot of the c alpha coordinates in the plane z=0
plt.scatter(calpha_atoms['X'], calpha_atoms['Y'])
plt.show()

#Tranform df into first 2 PCA components
pca = PCA(n_components=2)
pca_calpha_atoms = pd.DataFrame( data = pca.fit_transform(calpha_atoms),
                                columns = ['PCA1','PCA2'])

#Plot of the PCA components
plt.scatter(pca_calpha_atoms['PCA1'],pca_calpha_atoms['PCA2'])
plt.show()


#Rotate and mirror last plot to resemble 3D model
plt.scatter(-pca_calpha_atoms['PCA2'],-pca_calpha_atoms['PCA1'])
plt.xlim(-80,80)
plt.ylim(-90,70)
plt.show()

    

