import numpy as np
import dolfin
from dolfin import *

from mpi4py import MPI as pyMPI

from leopart import StokesStaticCondensation, FormsStokes
import geopart.stokes.incompressible
import dolfin_dg as dg

comm = pyMPI.COMM_WORLD
mpi_comm = MPI.comm_world


#load mesh,boundaries and coefficients from file
mark = {"Internal":0, "wall": 1,"inlet": 2,"outlet": 3 }

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

#output viscosity to paraview
XDMFFile(mpi_comm, "coeff_preview.xdmf").write_checkpoint(mu, "coeffs", 0)

#Define HDG element and function space
element_cls = geopart.stokes.incompressible.HDG2()
W = element_cls.function_space(mesh)

ds = dolfin.Measure('ds',domain=mesh,subdomain_data=boundaries)
n = dolfin.FacetNormal(mesh)

#Define boundary condition
U = element_cls.create_solution_variable(W)

p_in = dolfin.Constant(1.0)         # pressure inlet
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

#Form and solve Stokes
A, b = dolfin.PETScMatrix(), dolfin.PETScVector()
element_cls.solve_stokes(W, U, (A, b), weak_bcs, model)
uh, ph = element_cls.get_velocity(U), element_cls.get_pressure(U)

#Output solution p,u to paraview
dolfin.XDMFFile("pressure.xdmf").write_checkpoint(ph, "p")
dolfin.XDMFFile("velocity.xdmf").write_checkpoint(uh, "u")

flux = [dolfin.assemble(dolfin.dot(uh, n)*ds(i)) for i in range(len(mark))]

if comm.Get_rank() == 0:
    for key, value in mark.items():
        print("Flux_%s= %.15lf"%(key,flux[value]))