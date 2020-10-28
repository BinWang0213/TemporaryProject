import numpy as np
import geopart.stokes
import dolfin
import dolfin_dg as dg
import ufl
from ufl import sym, grad, Identity, div, dx
import matplotlib.pyplot as plt

#Load mesh and boundary
mesh = dolfin.Mesh()
hdf = dolfin.HDF5File(mesh.mpi_comm(), "mesh_boundaries.h5", "r")
hdf.read(mesh, "/mesh", False)

boundaries = dolfin.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)
mark = {"Internal":0, "wall": 1,"inlet": 2,"outlet": 3 }
hdf.read(boundaries, "/boundaries")
hdf.close()

#Define HDG element and function space
element_cls = geopart.stokes.HDG2()
W = element_cls.function_space(mesh)

ds = dolfin.Measure('ds',domain=mesh,subdomain_data=boundaries)
n = dolfin.FacetNormal(mesh)

#Define boundary condition
U = element_cls.solution_variable(W)

mu = dolfin.Constant(0.001)         # dynamic viscosity
p_in = dolfin.Constant(0.1)         # pressure inlet
p_out = dolfin.Constant(0.0)        # pressure outlet
noslip = dolfin.Constant([0.0]*mesh.geometry().dim())  # no-slip wall

#Boundary conditions
gN1 = (- p_out*dolfin.Identity(mesh.geometry().dim())) * n
Neumann_outlet=dg.DGNeumannBC(ds(mark["outlet"]), gN1)
gN2 = (- p_in*dolfin.Identity(mesh.geometry().dim())) * n
Neumann_inlet=dg.DGNeumannBC(ds(mark["inlet"]), gN2)
Dirichlet_wall=dg.DGDirichletBC(ds(mark["wall"]), noslip)

weak_bcs = [Dirichlet_wall,Neumann_inlet,Neumann_outlet]

#Body force term
f = dolfin.Constant([0.0]*mesh.geometry().dim())
model=geopart.stokes.StokesModel(eta=mu,f=f)

#Form Stokes
A, b = dolfin.PETScMatrix(), dolfin.PETScVector()
element_cls.solve_stokes(W, U, (A, b), weak_bcs, model)
uh, ph = element_cls.get_velocity(U), element_cls.get_pressure(U)


#visulize solution in Paraview
ufile_pvd = dolfin.File("velocity.pvd")
ufile_pvd << uh
pfile_pvd = dolfin.File("pressure.pvd")
pfile_pvd << ph #pressure is flipped for symmetric matrice

bd_pvd = dolfin.File("boundaries_mesh.pvd")
bd_pvd << boundaries


#Check boundary wall velocity

#Plot velocity u along bottom edge (0,0)->(12,0)
x=np.linspace(0,12,10)
vel_interp=np.array([uh([i,0.0]) for i in x])
print('BoundaryVelocity @ Line (0,0)->(12,0)=',np.sum(np.abs(vel_interp)))

plt.plot(x,vel_interp[:,0],label='Velocity x')
plt.plot(x,vel_interp[:,1],label='Velocity y')
plt.legend()
plt.show()