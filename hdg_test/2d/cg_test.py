import numpy as np
import dolfin
from dolfin import *
from mpi4py import MPI as pyMPI

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

#Define Taylor-Hood element and function space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

ds = dolfin.Measure('ds',domain=mesh,subdomain_data=boundaries)
n = dolfin.FacetNormal(mesh)

#Define boundary condition
p_in = dolfin.Constant(1.0)         # pressure inlet
p_out = dolfin.Constant(0.0)        # pressure outlet
noslip = dolfin.Constant([0.0]*mesh.geometry().dim())  # no-slip wall

#Boundary conditions

# No-slip Dirichlet boundary condition for velocity
bc0 = DirichletBC(W.sub(0), noslip, boundaries, mark["wall"])
bcs = [bc0]

#Neumann BC
gNeumann = - p_in  * inner(n, v) * ds(mark["inlet"]) + \
           - p_out * inner(n, v) * ds(mark["outlet"])

#Body force
f = Constant([0.0]*mesh.geometry().dim())

a = mu*inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx # The sign of the pressure has been flipped for symmetric system
L= inner(f, v)*dx + gNeumann

U = Function(W)
solve(a == L, U, bcs)

uh, ph = U.split()

#Output solution p,u to paraview
dolfin.XDMFFile("pressure.xdmf").write_checkpoint(ph, "p")
dolfin.XDMFFile("velocity.xdmf").write_checkpoint(uh, "u")

flux = [dolfin.assemble(dolfin.dot(uh, n)*ds(i)) for i in range(len(mark))]

if comm.Get_rank() == 0:
    for key, value in mark.items():
        print("Flux_%s= %.15lf"%(key,flux[value]))