import numpy as np
from BATL_generic import BATL_generic
class BATL_idl(BATL_generic):

	floatFmt = 'f' # 'f'->single, 'd'->double
	loadAll  = True
	varsToLoad = None

	def __init__(self, file, vars = None, precision= None):

		super(BATL_idl,self).__init__(file)

		if(vars is None):
			self.loadAll = True
		else:
			self.loadAll = False
			self.varsToLoad = vars

		if(precision is None):
			self.floatFmt = 'f'
		else:
			if(precision is 'single'):
				self.floatFmt = 'f'
			elif(precision is 'double'):
				self.floatFmt = 'd'
			else:
				print('Unknown precision option: {:s}', precision)
				return

		self.readBinary()

	#===============================================
	# Store data
	#===============================================
	def readBinary(self):
		import struct

		# Helper for processing Fortran record header/footers
		def processRecordTags(f,hdrsz):
			f.read(hdrsz)

		#===============================================
		#                   BEGIN
		#===============================================
		print('Loading binary file {:s}'.format(self.filename))

		# Size of record header/footer in bytes
		hdrsz    = 4 # bytes -> 32 bits
		charsz   = 1 # byte
		intsz    = 4
		singlesz = 4
		doublesz = 8
		floatsz  = singlesz if (self.floatFmt=='f') else doublesz
		npfltsz  = np.float32 if (self.floatFmt=='f') else np.float64

		#------------------------------------------------------------------------
		# Open file and perform read.
		#
		# This involves moving through the binary IDL output and processing it
		# appropriately. Note that Fortran binary files produce discrete
		# "records". These records are pre/appended with a header (and
		# identical footer) of 4 bytes which, stores the number of bytes in the
		# record (exlc header/footer). As a results, these tags must be accounted
		# for when advancing the file position.
		#------------------------------------------------------------------------
		with open(self.filename, 'rb') as f:

			# Read file header data
			processRecordTags(f,hdrsz)
			fileheader = struct.unpack(500*'c',f.read(charsz*500))
			processRecordTags(f,hdrsz)

			# Process step/time/dim/etc
			processRecordTags(f,hdrsz)
			self.step = struct.unpack('i',f.read(intsz))[0]
			self.time = struct.unpack(self.floatFmt,f.read(floatsz))[0]
			self.ndim = abs(struct.unpack('i',f.read(intsz))[0])
			self.npar = struct.unpack('i',f.read(intsz))[0]
			self.nvar = struct.unpack('i',f.read(intsz))[0]
			processRecordTags(f,hdrsz)

			print('step: {:d}'.format(self.step))
			print('time: {:e}'.format(self.time))
			print('ndim: {:d}'.format(self.ndim))
			print('npar: {:d}'.format(self.npar))
			print('nvar: {:d}'.format(self.nvar))

			# Process "zones/dim"
			processRecordTags(f,hdrsz)
			n_D = struct.unpack(3*'i',f.read(3*intsz))
			processRecordTags(f,hdrsz)

			# Read parameters
			if(self.npar>0):
				processRecordTags(f,hdrsz)
				parv = struct.unpack('f',f.read(singlesz))[0]
				processRecordTags(f,hdrsz)

			# Read variable/parameter names and remove
			# 'x','y','z' if necessary
			processRecordTags(f,hdrsz)
			varnames = (b''.join(struct.unpack_from(500*'s',f.read(500*charsz)))).decode('ascii')
			processRecordTags(f,hdrsz)

			varnames = varnames.strip(' ').split(' ')
			parnames = varnames[-self.npar:]
			varnames = varnames[abs(self.ndim):-self.npar]
			loadvars = varnames

			# Load grid points
			nZones = 1
			for i in range(0,abs(self.ndim)):
				nZones = nZones*n_D[i]
			processRecordTags(f,hdrsz)
			self.x = np.fromfile(f,dtype=np.dtype(np.float32),count=nZones)
			print('x bounds: {:e} {:e}'.format(np.amin(self.x),np.amax(self.x)))
			if(abs(self.ndim)>1):
				self.y = np.fromfile(f,dtype=np.dtype(np.float32),count=nZones)
				print('y bounds: {:e} {:e}'.format(np.amin(self.y),np.amax(self.y)))
			if(abs(self.ndim)>2):
				self.z = np.fromfile(f,dtype=np.dtype(np.float32),count=nZones)
				print('z bounds: {:e} {:e}'.format(np.amin(self.z),np.amax(self.z)))
			processRecordTags(f,hdrsz)

			# Determine which variables to read
			print('file variables: ', varnames)
			if(not self.loadAll):
				loadvars = self.varsToLoad
				print('variables to load: ', loadvars)

			# Read variable data
			for v in varnames:
				processRecordTags(f,hdrsz)
				if(v in loadvars):
					print('Loading {:s}'.format(v))
					self.data[v] = np.fromfile(f,dtype=np.dtype(npfltsz),count=nZones)
				else:
					f.seek(nZones*floatsz,1)
				processRecordTags(f,hdrsz)

			for v in loadvars:
				print('{:s} \n\tmin, max: {:e}  {:e}'.format(v,np.amin(self.data[v]),np.amax(self.data[v])))
