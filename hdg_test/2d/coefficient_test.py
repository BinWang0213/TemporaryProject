import numpy as np
import dolfin
from dolfin import *

#read mesh and boundaries from file
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "mesh_boundaries.h5", "r")
hdf.read(mesh, "/mesh", False)
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")
hdf.close()

#read viscosity coefficient from file
mu_ele = FunctionSpace(mesh, "DG", 0)
mu = Function(mu_ele)
hdf = HDF5File(mesh.mpi_comm(), "mesh_coeffs.h5", "r")
hdf.read(mu, "/mu")
hdf.close()

#Visulization !!!only for 2D
import matplotlib.tri as mtri
import matplotlib.pyplot as plt
plt.figure(figsize=(10,4))

nodes=mesh.coordinates()
cells=mesh.cells()

#plot mesh with coefficient
triang = mtri.Triangulation(nodes[:,0], nodes[:,1], triangles=cells)
plt.tripcolor(triang, facecolors=mu.vector()[:],linewidth=0.5,edgecolors='k',cmap='rainbow')
plt.colorbar()

plt.axis('equal')
plt.savefig('coefficient.png', bbox_inches='tight',dpi=100)