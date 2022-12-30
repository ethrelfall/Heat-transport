# this is example from https://www.firedrakeproject.org/demos/rayleigh-benard.py
# with small changes by Ed Threlfall, 30 December 2022


# Rayleigh-Benard Convection
# ==========================
# This problem involves a variable-temperature incompressible fluid.
# Variations in the fluid temperature are assumed to affect the momentum
# balance through a buoyant term (the Boussinesq approximation), leading
# to a Navier-Stokes equation with a nonlinear coupling to a
# convection-diffusion equation for temperature.
#
# We will set up the problem using Taylor-Hood elements for
# the Navier-Stokes part, and piecewise linear elements for the
# temperature. ::

from firedrake import *

# this boundary-refined mesh seems to give nicer results than the built-in uniform mesh
# mesh coords are warped in x and y according to 
# x -> (1/2)(1+atan(w*(x-1/2))/atan(w/2))

M = Mesh("refined_unit_square.msh")  # 40*40, stretched in x and y warp factor w=12

V = VectorFunctionSpace(M, "CG", 5)  # CG2, CG1 for pressure is Taylor-Hood
W = FunctionSpace(M, "CG", 4)
Q = FunctionSpace(M, "CG", 4)
Z = V * W * Q

upT = Function(Z)
u, p, T = split(upT)
v, q, S = TestFunctions(Z)

# Two key physical parameters are the Rayleigh number (Ra), which
# measures the ratio of energy from buoyant forces to viscous
# dissipation and heat conduction and the
# Prandtl number (Pr), which measures the ratio of viscosity to heat
# conduction. ::

x, y = SpatialCoordinate(M)

Ra = Constant(2e6)  # I have not been able to get much higher Ra values e.g. 3e6 to converge
Pr = Constant(10.0)

# Along with gravity, which points down. ::

g = Constant((0, 1))  # Ed: changed sign here because the buoyancy force is UP

F = (
    inner(grad(u), grad(v))*dx
    + inner(dot(grad(u), u), v)*dx
    - inner(p, div(v))*dx
    - (Ra/Pr)*inner(T*g, v)*dx
    + inner(div(u), q)*dx
    + inner(dot(grad(T), u), S)*dx
    + 1/Pr * inner(grad(T), grad(S))*dx
)

# There are two common versions of this problem.  In one case, heat is
# applied from bottom to top so that the temperature gradient is
# enforced parallel to the gravitation.  In this case, the temperature
# difference is applied horizontally, perpendicular to gravity.  It
# tends to make prettier pictures for low Rayleigh numbers, but also
# tends to take more Newton iterations since the coupling terms in the
# Jacobian are a bit stronger.  Switching to the first case would be a
# simple change of bits of the boundary associated with the second and
# third boundary conditions below::

bcs = [
    DirichletBC(Z.sub(0), Constant((0, 0)), (11, 12, 13, 14)),
    DirichletBC(Z.sub(2), Constant(1.0), (14,)),
    DirichletBC(Z.sub(2), Constant(0.0), (12,))
]

# Like Navier-Stokes, the pressure is only defined up to a constant.::

nullspace = MixedVectorSpaceBasis(
    Z, [Z.sub(0), VectorSpaceBasis(constant=True), Z.sub(2)])


# First off, we'll solve the full system using a direct solver.  As
# previously, we use MUMPS, so wrap the solve in ``try/except`` to avoid
# errors if it is not available. ::

from firedrake.petsc import PETSc

try:
   solve(F == 0, upT, bcs=bcs, nullspace=nullspace,
         solver_parameters={"mat_type": "aij",
                            "snes_monitor": None,
                            #"snes_view": None,
                            #"ksp_monitor_true_residual": None,
                            #"snes_converged_reason": None,
                            #"ksp_converged_reason": None,
                            "ksp_type": "gmres",
                            #"ksp_gmres_restart": 100,  #TRIALCODE
                            "pc_type": "lu",
                            "pc_factor_mat_solver_type": "mumps"})

except PETSc.Error as e:
    if e.ierr == 92:
        warning("MUMPS not installed, skipping direct solve")
    else:
        raise e


# do output here
u, p, T = upT.split()
u.rename("Velocity")
p.rename("Pressure")
T.rename("Temperature")
File("benard_mod.pvd").write(u, p , T)

# simple routine to evaluate Nusselt numbers (integral of heat flux at vertical walls)

normL = Function(V)
normL = Constant((-1.0,0.0))
fluxL = assemble(inner(normL, grad(T))*ds(14))
print('flux on LHS is:')
print(fluxL)
normR = Function(V)
normR = Constant((1.0,0.0))
fluxR = assemble(inner(normR, grad(T))*ds(12))
print('flux on RHS is:')
print(fluxR)

# TRIALCODE run again - commented out as this crashes immediately (I'm sure it's my fault)
# why can't I just run the thing again with the previous solution as initial guess?
# do I need to instantiate something and then re-use it, instead of calling solve twice?

#Ra = Constant(3e6)  # idea is to set a larger Ra but use output from previous as initial guess

#F = (
#    inner(grad(u), grad(v))*dx
#    + inner(dot(grad(u), u), v)*dx
#    - inner(p, div(v))*dx
#    - (Ra/Pr)*inner(T*g, v)*dx
#    + inner(div(u), q)*dx
#    + inner(dot(grad(T), u), S)*dx
#    + 1/Pr * inner(grad(T), grad(S))*dx
#)

#solve(F == 0, upT, bcs=bcs, nullspace=nullspace,
#      solver_parameters={"mat_type": "aij",
#                         "snes_monitor": None,
#                         #"snes_view": None,
#                         #"ksp_monitor_true_residual": None,
#                         #"snes_converged_reason": None,
#                         #"ksp_converged_reason": None,
#                         "ksp_type": "gmres",
#                         #"ksp_gmres_restart": 100,  #TRIALCODE
#                         "pc_type": "lu",
#                         "pc_factor_mat_solver_type": "mumps"})

# outputs ...
#u, p, T = upT.split()
#u.rename("Velocity")
#p.rename("Pressure")
#T.rename("Temperature")
#File("benard_mod2.pvd").write(u, p , T)

print('finished - quitting.')
quit()

# I've deleted the rest of the example as trying to understand why problem
# does not converge for exact solver when Ra > approx 2e6


