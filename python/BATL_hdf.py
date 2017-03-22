import h5py
import numpy as np
from BATL_generic import BATL_generic
class BATL_hdf(BATL_generic):

	loadAll  = True
	varsToLoad = None

	def __init__(self, file, vars = None):

		super(BATL_hdf,self).__init__(file)

		if(vars is None):
			self.loadAll = True
		else:
			self.loadAll = False
			self.varsToLoad = vars

		self.readHDF5()

	#===============================================
	# Store data
	#===============================================
	def readHDF5(self):
		import struct

		#===============================================
		#                   BEGIN
		#===============================================
		print('Loading HDF5 file {:s}'.format(self.filename))
		with h5py.File(self.filename,"r") as f:



			print('var names: ', f['plotVarNames'].value)

			blkcenter = f['coordinates'].value
			blkbnds   = f['bounding box'].value
			nI        = f['Integer Sim Metadata'].value[0]
			nJ        = f['Integer Sim Metadata'].value[1]
			nK        = f['Integer Sim Metadata'].value[2]
			print(blkcenter.shape)
			print(blkbnds.shape)
			print(nI, nJ, nK)


			# nb,nx,ny,nz = f['posx'].shape
			# self.x = f['posx'].reshape(nb*nx*ny*nz)
			# self.y = f['posy'].reshape(nb*nx*ny*nz)
			# self.z = f['posz'].reshape(nb*nx*ny*nz)