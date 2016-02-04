"""
Solve (u_xx + u_yy)*nu_u = u_t on a uniform grid
and (v_xx + v_yy)*nu_v = v_t on a uniform grid

with known solution

u, v = exp(-nu*pi^2*(a^2+b^2)t)*sin(a*pi*x)*sin(b*pi*y)
and initial condition
u, v = sin(a*pi*x)*sin(b*pi*y)
for integer choices a and b
this allows for a lot of flexibility in how things want to be defined

lets pick u_a = 1, u_b = 1, nu_u = 1
and v_a = 4, v_b = 2, nu_v = 0.5
and tf = 1
python diffusion2D.py n
for the number of cells in each dimension

Run this problem as python poisson2D.py N
where N is the number of cells. The code does a rough accuracy estimation,
and outputs the accuracy for both u and v. The accuracy should be around 2 if 
the code is performing correctly.

"""
import math as math
import sys
import src.base.blocks as b
import src.base.flux as f
import src.base.problem as p
import src.base.source as s
import numpy as np

""" 
lets define a uniform square mesh on [-1, 1] x [1, 1]
and create boundary blocks as we go,
initializing based on the exact solution, and naming the block by its coordinates
"""
u_a = 1
u_b = 1
nu_u = 1
v_a = 1
v_b = 1
nu_v = 0.5
tf = 1

""" Sources defined here """
def time(B,P):
	u = math.exp(-nu_u*math.pi*math.pi*(u_a*u_a+u_b*u_b)*B.t)*math.sin(u_a*math.pi*P['x'])*math.sin(u_b*math.pi*P['y'])
	v = math.exp(-nu_v*math.pi*math.pi*(v_a*v_a+v_b*v_b)*B.t)*math.sin(v_a*math.pi*P['x'])*math.sin(v_b*math.pi*P['y'])
	return {'u':u,'v':v}

""" Fluxes defined here """
def difference(B,N,G):
	return dict((s,(N[s]-B[s])/G['d']) for s in B.state)

def diff2D(N):
	d = 2./float(N) # spacing, delta 
	# initialize with exact solution at t = 0
	B = [b.Block('('+str(i*d-d/2-1)+','+str(j*d-d/2-1)+')',
		{'u':math.sin(u_a*math.pi*(i*d-d/2-1))*math.sin(u_b*math.pi*(j*d-d/2-1)),
		 'v':math.sin(v_a*math.pi*(i*d-d/2-1))*math.sin(v_b*math.pi*(j*d-d/2-1))}) \
		for i in range(0,N+2) for j in range(0,N+2)]

	# Flux geometry	
	G = {'type':'edge','d':d*d,'m':[]}
	n = N+2 # add two for the boundaries
	for i in range(1,n-1):
		for j in range(1,n-1):
			# Add fluxes, no "normals" so we cheat and define them cleverly
			for k in [(i-1)*n+j, i*n+j-1,(i+1)*n+j, i*n+j+1]:
				B[i*n+j].addFlux(f.Flux(B[k],difference,G))
			B[i*n+j]['u'] = 0.
			B[i*n+j]['v'] = 0.

	interiorBlocks = [B[i*n+j] for i in range(1,n-1) for j in range(1,n-1)]
	boundaryBlocks = []
	bcRange = [j         for j in range(1,n-1)] + \
						[(n-1)*n+j for j in range(1,n-1)] + \
						[i*n       for i in range(0,n)] + \
						[i*n+n-1   for i in range(0,n)] 

	for k in bcRange:
		(x,y) = eval(B[k].name)
		B[k].addSource(s.Source(time,{'x':x,'y':y},B[k].name))
		boundaryBlocks.append(B[k])

	# solve the problem on the interior blocks
	P = p.Problem(interiorBlocks,boundaryBlocks)
	P.solveUnst(np.linspace(0,tf,10))
	# P.printSolution()
	Eu = 0
	Ev = 0
	t = tf
	for bb in interiorBlocks:
		(x,y) = eval(bb.name)
		ue = math.exp(-nu_u*math.pi*math.pi*(u_a*u_a+u_b*u_b)*t)*math.sin(u_a*math.pi*x)*math.sin(u_b*math.pi*y)
		ve = math.exp(-nu_v*math.pi*math.pi*(v_a*v_a+v_b*v_b)*t)*math.sin(v_a*math.pi*x)*math.sin(v_b*math.pi*y)
		Eu += (bb['u']-ue)**2
		Ev += (bb['v']-ve)**2

	# Eu = math.sqrt(sum([(math.exp(-nu_u*math.pi*math.pi*(u_a*u_a+u_b*u_b)*tf)*math.sin(u_a*math.pi*block.name[0])*math.sin(u_b*math.pi*block.name[1]) \
	# 	-block.state['u'])**2 for block in interiorBlocks])/(n-2)/(n-2))
	# Ev = math.sqrt(sum([(math.exp(-nu_v*math.pi*math.pi*(v_a*v_a+v_b*v_b)*tf)*math.sin(v_a*math.pi*block.name[0])*math.sin(v_b*math.pi*block.name[1]) \
		# -block.state['v'])**2 for block in interiorBlocks])/(n-2)/(n-2))
	return (math.sqrt(Eu/(n-2)/(n-2)),math.sqrt(Ev)/(n-2)/(n-2))

def test():
	n = 3
	Error = [diff2D(n),diff2D(n*2)]
	Rate = [(math.log(Error[1][0])-math.log(Error[0][0]))/(math.log(2./(2*n))-math.log(2./(n))),
	(math.log(Error[1][1])-math.log(Error[0][1]))/(math.log(2./(2*n))-math.log(2./(n)))]
	return Rate
