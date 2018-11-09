#!/usr/bin/env python
# coding: utf-8
"""
Created Mon Nov 5 20:16:36 EST 2018
@author: ShengMao

This file uses modified Riks method to solve the wrinkling problem:
solve: {(u, g), F(u, g) == 0},
"""
##################################################################################
# import pkgs
##################################################################################
# standard pkgs from the envs
import dolfin as dl
import matplotlib.pyplot as plt
import numpy as np
import sys
# in-house code
from mesh_quad import gen_mesh
# parameter settings
deg = 4
# set the number of integration points
dl.parameters["form_compiler"]["quadrature_degree"] = deg
# set the font size
plt.rcParams.update({'font.size': 15})
# tolerance
tol = 1e-8
##################################################################################
# set the geometry accordingly
##################################################################################
"""
define the class of geometry of a rectangular with length height and
a film on top
"""
class Geometry():
    def __init__(self, length, height, filmHeight):
        # length, total thickness,
        self.length = length
        self.height = height
        self.filmHeight = filmHeight
    
    def gen_mesh(self, numElementsFilm):
        gen_mesh(self.length, self.height, self.filmHeight, numElementsFilm)

    def read_mesh(self, numElementsFilm, reGen=True):
        if reGen:
            self.gen_mesh(numElementsFilm=2.0)
        fileName = "mesh/geometry_quad.xdmf"
        self.mesh = dl.Mesh()
        f = dl.XDMFFile(dl.mpi_comm_world(), fileName)
        f.read(self.mesh)
        f.close()


# Sub domain for Periodic boundary condition
class PeriodicBoundary(dl.SubDomain):
    """
    create periodic boundary with domain info, only left to right
    """
    def __init__(self, domain):
        dl.SubDomain.__init__(self)
        self.L = domain.length

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return dl.near(x[0], 0) and on_boundary

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - self.L
        y[1] = x[1]


##################################################################################
# set the physical parameters
##################################################################################
"""
define the class containing all the physical parameters:
    -- Lame constants of the substrate: muSub, lmbdaSub
    -- modulus ration: ratFilmSub
"""
class Physical_Params():
    def __init__(self, lmbdaSub, muSub, ratFilmSub):
        # length, total thickness,
        self.lmbdaSub = lmbdaSub
        self.muSub = muSub
        self.ratFilmSub = ratFilmSub
    
    # create domain dependent constants
    def create_Kfilm_Ksub(self, ratFilmSub, Ksub, domain):
        Ht, Hf = domain.height, domain.filmHeight
        return dl.Expression('x[1] >= Ht - Hf - tol ? k*Ks : Ks', degree=0,
                             tol = tol, k = ratFilmSub, Ks = Ksub, Hf = Hf, Ht = Ht)
    
    # fill in lame constants
    def create_Lame(self, geometry):
        ratFilmSub, lmbdaSub, muSub = self.ratFilmSub, self.lmbdaSub, self.muSub
        self.lmbda = self.create_Kfilm_Ksub(ratFilmSub, lmbdaSub, geometry)
        self.mu = self.create_Kfilm_Ksub(ratFilmSub, muSub, geometry)

##################################################################################
# create specific elements
##################################################################################
def create_mixed_element(domain, order=1):
    VecElement = dl.VectorElement("CG", domain.mesh.ufl_cell(), order)
    FinElement = dl.FiniteElement("CG", domain.mesh.ufl_cell(), order)
    return dl.MixedElement([FinElement, VecElement])


def create_vector_element(domain, order=1):
    VecElement = dl.VectorElement("CG", domain.mesh.ufl_cell(), order)
    return VecElement

##################################################################################
# growth factor
##################################################################################
# this is repeating the create_Kfilm_Ksub function! need to rethink about it!
def create_growthFactor(domain, filmGrowth, subGrowth, growthInitial=0.0):
    """
        change gF to change the growth rate
    """
    # obtain the geometry parameters
    length, height, filmHeight = domain.length, domain.height, domain.filmHeight
    kFilm = filmGrowth; kSub = subGrowth
    growthFactor = dl.Expression('x[1] >= H - Hf - tol ? kf*gF : ks*gF',
                                 degree=0, tol=tol, gF=growthInitial,
                                 kf=filmGrowth, ks=subGrowth,
                                 Hf=filmHeight, H=height)
    return growthFactor


# translate the growth rate to matrix form
# currently only allowing uniaxial/isotropic growth
def create_Fg(growthFactor, type="uniaxial"):
    gf = growthFactor
    if type == "uniaxial":
        Fg = dl.as_matrix(((1.0 + gf, 0.0), (0., 1.0)))
    elif type == "isotropic":
        Fg = dl.as_matrix(((1.0 + gf, 0.0), (0., 1.0 + gf)))
    else:
        Fg = None
        print ("Error specifying growth type")
        sys.exit()
    return Fg

##################################################################################
# basic kinematics
##################################################################################
def cal_neoHookean(F, physParams):
    """
        neoHookean elastic energy density: psi
    """
    lmbda, mu = physParams.lmbda, physParams.mu
    ln = dl.ln
    J  = dl.det(F)
    Ic = dl.tr(F.T*F)     # Invariants of deformation tensors
    return (lmbda/4)*(J**2 - 2*ln(J) - 1) + (mu/2)*(Ic-2) - mu*ln(J)


class BucklingProblem(dl.OptimisationProblem):
    """
        Optimization problem for wrinkling analysis
    """
    def __init__(self, w, Energy, Residual, Jacobian):
        dl.OptimisationProblem.__init__(self)
        self.w = w
        self.Energy = Energy
        self.Residual = Residual
        self.Jacobian = Jacobian
    # Objective function
    def f(self, x):
        self.w.vector()[:] = x
        return dl.assemble(self.Energy)
    # Gradient of the objective function
    def F(self, b, x):
        self.w.vector()[:] = x
        dl.assemble(self.Residual, tensor=b)
    # Hessian of the objective function
    def J(self, A, x):
        self.w.vector()[:] = x
        dl.assemble(self.Jacobian, tensor=A)

# create upper bounds and lower bounds of displacement
def create_bounds(wMin, wMax, functionSpace, boundaryCondtions=None):
    upperBound = dl.Expression(("w0max", "w1max"), w0max=wMax[0], w1max=wMax[1], degree=0)
    lowerBound = dl.Expression(("w0min", "w1min"), w0min=wMin[0], w1min=wMin[1], degree=0)
    wUpperBound = dl.interpolate(upperBound, functionSpace)
    wLowerBound = dl.interpolate(lowerBound, functionSpace)
    
    # apply boundary conditions to the bounds
    if boundaryCondtions:
        for bc in boundaryCondtions:
            bc.apply(wUpperBound.vector())
            bc.apply(wLowerBound.vector())
    return wLowerBound, wUpperBound

##################################################################################
# main function
##################################################################################
def main(dt, tEnd,
         length=10.0, height=5.0, numElementsFilm=1.5, reGen=True,
         ratLame=5.0, ratFilmSub=100.0):
    fileDir = ("results-dt-%.2f-tEnd-%.0f-L-%.1f-H-%.1f-ratioLame-%.0f-ratioFilm-%.0f" %
               (dt, tEnd, length, height, ratLame, ratFilmSub))
    # create geometry
    rectDomain = Geometry(length=length, height=height, filmHeight=0.05)
    rectDomain.read_mesh(numElementsFilm, reGen)
    # define periodic boundary conditions
    periodicBC = PeriodicBoundary(rectDomain)
    # specify physical parameters
    physParams = Physical_Params(lmbdaSub=ratLame, muSub=1.0, ratFilmSub=100.0)
    physParams.create_Lame(rectDomain)
    # create discrete function space
    element = create_vector_element(rectDomain, order=1)
    W = dl.FunctionSpace(rectDomain.mesh, element, constrained_domain=periodicBC)
    
    # define boundaries
    def left(x, on_boundary):
        return dl.near(x[0], 0.) and on_boundary

    def right(x, on_boundary):
        return dl.near(x[0], rectDomain.length) and on_boundary

    def bottom(x, on_boundary):
        return dl.near(x[1], 0.0) and on_boundary

    def top(x, on_boundary):
        return dl.near(x[1], rectDomain.height) and on_boundary

    def corner(x, on_boundary):
        return dl.near(x[0], 0.0) and dl.near(x[1], 0.0)
    
    # define fixed boundary
    bcs = [dl.DirichletBC(W.sub(0), dl.Constant(0), left),
           dl.DirichletBC(W.sub(0), dl.Constant(0), right),
           dl.DirichletBC(W.sub(1), dl.Constant(0), bottom),
           dl.DirichletBC(W.sub(0), dl.Constant(0), corner, method="pointwise")]
    bcs=bcs[-2:]
    # the variable to solve for
    w = dl.Function(W, name="Variables at current step")
    # test and trial function
    dw = dl.TrialFunction(W); w_ = dl.TestFunction(W)
    # dealing with physics of growth
    growthFactor = create_growthFactor(rectDomain, filmGrowth=1, subGrowth=0)
    Fg = create_Fg(growthFactor, "uniaxial")
    # kinematics
    I = dl.Identity(2); F = I + dl.grad(w); Fe = F * dl.inv(Fg)
    # write the variational form from potential energy
    psi = cal_neoHookean(Fe, physParams)
    Energy = psi*dl.dx
    Residual = dl.derivative(Energy, w, w_)
    Jacobian = dl.derivative(Residual, w, dw)
    problem = BucklingProblem(w, Energy, Residual, Jacobian)
    wLowerBound, wUpperBound = create_bounds(wMin = [-0.5, -0.5],
                                             wMax = [0.5, 0.5],
                                             functionSpace = W,
                                             boundaryCondtions = bcs)
    # Create the PETScTAOSolver
    solver = dl.PETScTAOSolver()
    TAOSolverParameters = {"method": "tron", "maximum_iterations":1000,
                           "monitor_convergence": True}
    solver.parameters.update(TAOSolverParameters)
    #solver.parameters["report"] = False
    #solver.parameters["linear_solver"] = "umfpack"
    #solver.parameters["line_search"] = "gpcg"
    #solver.parameters["preconditioner"] = "ml_amg"

    growthRate= 0.1;growthEnd = growthRate * tEnd
    gF = 0.0
    outDisp = dl.File(fileDir + "/displacement.pvd")
    outDisp << (w, 0.0)

    while gF <= growthEnd - tol:
        gF += growthRate * dt
        growthFactor.gF = gF
        nIters, converged = solver.solve(problem, w.vector(),
                                         wLowerBound.vector(),
                                         wUpperBound.vector())
        w.vector()[:] += 2e-04 * np.random.uniform(-1, 1, w.vector().local_size())
        print("--- growth = %.4f, niters = %d, dt = %.5f -----" % (1+gF, nIters, dt))
        outDisp << (w, gF)

################################ test ########################################
"""
V = dl.FunctionSpace(rectDomain.mesh, "CG", 1)
# test Lame constants
lmbda = dl.Function(V)
lmbda.interpolate(physParams.lmbda)
out = dl.File("test/lame.pvd")
out << (lmbda, 0.0)
# test growth factors
gF = dl.Function(V)
gF.interpolate(growthFactor)
out = dl.File("test/growth.pvd")
out << (gF, 0.0)
growthFactor.gF = 0.01

#w.interpolate(dl.Expression(('sin(2*pi*x[0]/L)', 'cos(2*pi*x[0]/L)'),
#                             L = rectDomain.length, pi=np.pi, degree=1))

print(dl.assemble(Energy))
"""


