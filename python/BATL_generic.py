import numpy as np
import scipy.interpolate
class BATL_generic(object):

	def __init__(self, file):
		self.filename = file
		#--------------
		# DATA
		#--------------
		self.x = None
		self.y = None
		self.z = None
		self.paramters = dict()
		self.data = dict()
		self.step = -1
		self.time = -1
		self.ndim = -1
		self.npar = -1
		self.nvar = -1

	def getUniform(self, var, bndBox=None, deltas=None, nPts=None, xGrid=None, yGrid=None, zGrid=None):
		from math import ceil

		# Construct bounding box
		bb = [[0.0,0.0],[0.0,0.0],[0.0,0.0]]
		if(bndBox is None):
			if(self.x is not None):
				bb[0][0] = np.amin(self.x)
				bb[0][1] = np.amax(self.x)
			if(self.y is not None):
				bb[1][0] = np.amin(self.y)
				bb[1][1] = np.amax(self.y)
			if(self.z is not None):
				bb[2][0] = np.amin(self.z)
				bb[2][1] = np.amax(self.z)
		else:
			bb = bndBox

		# Construct number of sample points in each dimension
		n_d = [1,1,1]
		if(deltas is not None):
			n_d[0] = math.ceil((bb[0][1]-bb[0][1])/float(deltas[0]))
			n_d[1] = math.ceil((bb[1][1]-bb[1][1])/float(deltas[1])) if(self.ndim>1) else 1
			n_d[2] = math.ceil((bb[2][1]-bb[2][1])/float(deltas[2])) if(self.ndim>2) else 1
		elif(nPts is not None):
			n_d = nPts
			for i in range(self.ndim,3):
				n_d[i]=1


		X = None
		Y = None
		Z = None
		if(self.ndim==1):
			if(xGrid is not None):
				if(np.all(np.diff(xGrid) > 0)):
					X = xGrid
				else:
					print('[getUniform] Provided x grid is not strictly increasing')
			else:
				X = np.linspace(bb[0][0],bb[0][1],n_d[0])

			intfunc = scipy.interpolate.interp1d(self.xs,self.data[var],'nearest')
			F = intfunc(X)
			return X,Y,Z,F

		elif(self.ndim==2):
			if(xGrid and yGrid):
				X = xGrid
				Y = xGrid
			else:
				xs  = np.linspace(bb[0][0],bb[0][1],n_d[0])
				ys  = np.linspace(bb[1][0],bb[1][1],n_d[1])
				X,Y = np.meshgrid(xs,ys,indexing='xy')

			F = scipy.interpolate.griddata((self.x,self.y),self.data[var],(X,Y),method='nearest',fill_value=-777.0)
			return X,Y,Z,F
		elif(self.ndim==3):
			print('[getUniform] 3D not implemented')
			return None, None, None, None
		else:
			print('[getUniform] Unhandled ndim value {:d}'.format(self.ndim))
			return None, None, None, None



















