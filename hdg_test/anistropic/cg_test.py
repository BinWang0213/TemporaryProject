import numpy as np
import dolfin
from dolfin import *
from mpi4py import MPI as pyMPI

comm = pyMPI.COMM_WORLD
mpi_comm = MPI.comm_world

#mark whole boundary, inflow and outflow will overwrite)
class Noslip(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1.0)

#Create a unit box mesh
n_ele = 6
aspect_ratio = 3

mesh = BoxMesh(comm, Point(0.0, 0.0,0.0), Point(1.0, 1.0, 1.0),
                     n_ele, n_ele, n_ele*aspect_ratio)

#read mesh and boundaries from file
boundaries = MeshFunction('size_t', mesh, mesh.topology().dim() - 1)

mark = {"Internal":0, "wall": 1,"inlet": 2,"outlet": 3 }

boundaries.set_all(mark["Internal"])
wall=Noslip()
wall.mark(boundaries, mark["wall"])
left = Left()
left.mark(boundaries, mark["inlet"])
right = Right()
right.mark(boundaries, mark["outlet"])

#read viscosity coefficient from file
mu = Constant(0.001)

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

#Weak form
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