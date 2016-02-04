"""
Solve u_xx + u_yy = f(x,y) on a uniform grid
and v_xx + v_yy = f(x,y) on a uniform grid

with known solution

u = exp(xy) 	 				f['u'] = (x^2+y^2)*exp(xy)
v = exp((x^2+y^2))  f['v'] = (4*(x^2+y^2+1))*exp((x^2+y^2)) 

This is discretized

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


""" Fluxes defined here """
def difference(B,N,P):
	return dict((s,(N[s]-B[s])/P['d']) for s in B.state)

def poisson2D(N):
	""" 
	lets define a uniform square mesh on [-1, 1] x [1, 1]
	and create boundary blocks as we go,
	initializing based on the exact solution, and naming the block by its coordinates
	"""
	d = 2./float(N) # spacing, delta X
	B = [b.Block('('+str(i*d-d/2-1)+','+str(j*d-d/2-1)+')',{'u':math.exp((i*d-d/2-1)*(j*d-d/2-1)),
		'v':math.exp((i*d-d/2-1)**2+(j*d-d/2-1)**2)}) for i in range(0,N+2) for j in range(0,N+2)]

	# interior sources
	# use the name to get the source values
	for block in B:
		(x,y) = eval(block.name)
		block.addSource(s.Source(s.constant,{'u':-(x*x+y*y)*math.exp(x*y),
			'v':-4.0*(x*x+y*y+1.0)*math.exp(x*x+y*y)},'constant'))
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

	# solve the problem on the interior blocks
	P = p.Problem(interiorBlocks)
	P.solve()
	Eu = math.sqrt(sum([(math.exp(eval(block.name)[0]*eval(block.name)[1])-block['u'])**2 for block in interiorBlocks])/(n-2)/(n-2))
	Ev = math.sqrt(sum([(math.exp(eval(block.name)[0]**2+eval(block.name)[1]**2)-block['v'])**2 for block in interiorBlocks])/(n-2)/(n-2))
	return (Eu,Ev)

def test():
	n = 3
	Error = [poisson2D(n),poisson2D(n*2)]
	Rate = [(math.log(Error[1][0])-math.log(Error[0][0]))/(math.log(2./(2*n))-math.log(2./(n))),
	(math.log(Error[1][1])-math.log(Error[0][1]))/(math.log(2./(2*n))-math.log(2./(n)))]
	return Rate
