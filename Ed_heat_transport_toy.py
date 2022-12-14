from firedrake import *

# proxy for "coupled" problem in conductive heat transport
# left-hand-side unit square is nominally beam region
# right-hand-side is nominally reactor wall, cooled to const temperatures

# this needs:
# 1) turning into two coupled problems (right now is tightly-coupled)
# 2) add radiation i.e. black body \sigma T**4 (Stefan's Law)

mesh=Mesh("bi_unit_square.msh")

V = FunctionSpace(mesh, "CG", 6)

u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(mesh)

# source function
# pick something without sharp edges (I used "Frankenstein's Monster"-type curve)
bf = 1+conditional(And(ge(y, -0.25), le(y, 0.25)), exp(16-1/(0.25*0.25-y**2)), 0)

bcLL = DirichletBC(V, bf, 16)  # source on LHS edge
bcR  = DirichletBC(V, 0.0, (18,19,20)) # const temperature

theta = 0.0*pi/180 # nonzero angle may cause artifacts for some BC choices
bhat=as_vector([cos(theta),sin(theta)])

k_par = conditional(le(x,0.5), 50.0, 1.0)
k_per = conditional(le(x,0.5), 0.0,  1.0)

flux = k_par * bhat * dot(bhat, grad(u)) + k_per * (grad(u) - bhat * dot(bhat, grad(u)))

a = inner(flux, grad(v))*dx
f = Function(V)
f.interpolate(0.0*x)
L = inner(f,v)*dx

T = Function(V)

solve( a==L, T, bcs=[bcLL, bcR])

# QOIs can be T in centre of domain and distribution of heat flux out (total is of course equal to heat flux in through LHS)
# T in centre basically given by T on left, T on right, and balance between plasma and metal heat conductivities
# cf. 1D model T_mid = T_L (k_L/(k_L+k_R)) + T_R (k_R/(k_L+k_R)) so if k_L >> k_R T_mid -> T_L
File("Ed_heat_transport_toy.pvd").write(T)  
