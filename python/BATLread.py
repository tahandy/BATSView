import h5py
import numpy as np
from scipy.interpolate import *
import matplotlib.pyplot as plt
import glob

def getData(inputFile,field):
	if field == 'flux':
		tr = inputFile['trkev'].value
		tr = 1000.0*tr/100.
		data = 2.36e23*tr*tr*tr
	elif field == 'front_velocity':
		tr = inputFile['trkev'].value
		tr = 1000.0*tr
		flux = 2.36e23*tr*tr*tr
		dens = inputFile['rho'].value
		dens = 2.0*dens / (14.0*1.67e-24)
		data = flux/dens
	elif field == "optical_depth":
		#Assume a cross section of 0.5*1.0e-18
		#I just need the density then
		data = inputFile['rho'].value / (14.*1.67e-24)
	elif field =='NI':
		zavg  = np.array(inputFile['zavg'].value)
		dens  = np.array(inputFile['rho'].value)/(14.*1.67e-24)
		ones  = np.ones( zavg.shape)
		zeros = np.zeros( zavg.shape)
		temp = np.where( (zavg>=0.0) & (zavg < 1.0),ones,zeros)
		#temp = np.where( zavg>=0.0, ones,zeros)
		#temp = np.where( zavg<1.0, temp,zeros)
		NI    = np.zeros( zavg.shape )
		NI    = (1.0-zavg)*temp#*dens
		return NI
	elif field[:-2] == "alpha":
		#First I need the relative populations for all stages of N
		#Need number densities.
		zavg  = np.array(inputFile['zavg'].value)
		#nb,nx,ny,nz = zavg.shape
		#zavg = zavg.reshape(nb*nx*ny*nz)
		dens  = np.array(inputFile['rho'].value)/(14.*1.67e-24)
		NI    = np.zeros( zavg.shape )
		NII   = np.zeros( zavg.shape )
		NIII  = np.zeros( zavg.shape )
		NIV   = np.zeros( zavg.shape )
		NV    = np.zeros( zavg.shape )
		NVI   = np.zeros( zavg.shape )
		NVII  = np.zeros( zavg.shape )
		NVIII = np.zeros( zavg.shape )
		ones  = np.ones( zavg.shape)
		zeros = np.zeros( zavg.shape)
		#Between zero and one
		temp = np.where( (zavg>=0.0) & (zavg < 1.0),ones,zeros)
		NI  = NI  + (1.0-zavg)*temp*dens
		NII = NII + (zavg-0.0)*temp*dens
		#Between one and two
		temp = np.where( (zavg>=1.0) & (zavg < 2.0),ones,zeros)
		NII  = NII  + (2.0-zavg)*temp*dens
		NIII = NIII + (zavg-1.0)*temp*dens
		#Between two and three
		temp = np.where( (zavg>=2.0) & (zavg < 3.0),ones,zeros)
		NIII  = NIII + (3.0-zavg)*temp*dens
		NIV   = NIV  + (zavg-2.0)*temp*dens
		#Between three and four
		temp = np.where( (zavg>=3.0) & (zavg < 4.0),ones,zeros)
		NIV  = NIV   + (4.0-zavg)*temp*dens
		NV   = NV    + (zavg-3.0)*temp*dens
		#Between four and five
		temp = np.where( (zavg>=4.0) & (zavg < 5.0),ones,zeros)
		NV    = NV   + (5.0-zavg)*temp*dens
		NVI   = NVI  + (zavg-4.0)*temp*dens
		#Between five and six
		temp = np.where( (zavg>=5.0) & (zavg < 6.0),ones,zeros)
		NVI    = NVI   + (6.0-zavg)*temp*dens
		NVII   = NVII  + (zavg-5.0)*temp*dens
		#Between six and seven
		temp = np.where( (zavg>=6.0) & (zavg < 7.0),ones,zeros)
		NVII   = NVII  + (7.0-zavg)*temp*dens
		NVIII  = NVIII + (zavg-6.0)*temp*dens

		#Recombination Rates
		T = np.array(inputFile['tikev'].value)*1000.0
		R21 = 6.387E-10/(np.sqrt(T/9.467E-02)*np.power(1.0+np.sqrt(T/9.467E-02),1.0-(0.7308+0.2440*np.exp(-6.739E+04/T)))*np.power(1.0+np.sqrt(T/2.954E+06),1.0+(0.7308+0.2440*np.exp(-6.739E+04/T)))) + np.power(T,-1.5)*(1.658E-08*np.exp(-1.265E+01/T)+2.760E-08*np.exp(-8.425E+01/T)+2.391E-09*np.exp(-2.964E+02/T)+7.585E-07*np.exp(-5.923E+03/T)+3.012E-04*np.exp(-1.278E+05/T)+7.132E-04*np.exp(-2.184E+05/T))
		R32 = 2.410E-09/(np.sqrt(T/1.231E-01)*np.power(1.0+np.sqrt(T/1.231E-01),1.0-(0.7948+0.0774*np.exp(-1.016E+05/T)))*np.power(1.0+np.sqrt(T/3.016E+06),1.0+(0.7948+0.0774*np.exp(-1.016E+05/T)))) + np.power(T,-1.5)*(7.712E-08*np.exp(-7.113E+01/T)+4.839E-08*np.exp(-2.765E+02/T)+2.218E-06*np.exp(-1.439E+04/T)+1.536E-03*np.exp(-1.347E+05/T)+3.647E-03*np.exp(-2.496E+05/T)+4.234E-05*np.exp(-2.204E+06/T))
		R43 = 7.923E-10/(np.sqrt(T/3.750E+00)*np.power(1.0+np.sqrt(T/3.750E+00),1.0-(0.7768+0.0223*np.exp(-7.206E+04/T)))*np.power(1.0+np.sqrt(T/3.468E+06),1.0+(0.7768+0.0223*np.exp(-7.206E+04/T)))) + np.power(T,-1.5)*(3.386E-06*np.exp(-1.406E+03/T)+3.036E-05*np.exp(-6.965E+03/T)+5.945E-05*np.exp(-2.604E+04/T)+1.195E-03*np.exp(-1.304E+05/T)+6.462E-03*np.exp(-1.965E+05/T)+1.358E-03*np.exp(-4.466E+06/T))
		R54 = 1.533E-10/(np.sqrt(T/1.823E+02)*np.power(1.0+np.sqrt(T/1.823E+02),1.0-0.6682)*np.power(1.0+np.sqrt(T/7.751E+06),1.0+0.6682)) + np.power(T,-1.5)*(2.040E-06*np.exp(-3.084E+03/T)+6.986E-05*np.exp(-1.332E+04/T)+3.168E-04*np.exp(-6.475E+04/T)+4.353E-03*np.exp(-1.181E+05/T)+7.765E-04*np.exp(-6.687E+05/T)+5.101E-03*np.exp(-4.778E+06/T))
		R65 = 6.245E-11/(np.sqrt(T/1.957E+03)*np.power(1.0+np.sqrt(T/1.957E+03),1.0-0.4985)*np.power(1.0+np.sqrt(T/2.177E+07),1.0+0.4985)) + np.power(T,-1.5)*(5.761E-03*np.exp(-3.860E+06/T)+3.434E-02*np.exp(-4.883E+06/T)-1.660E-03*np.exp(-6.259E+06/T))
		R76 = 2.388E-10/(np.sqrt(T/3.960E+02)*np.power(1.0+np.sqrt(T/3.960E+02),1.0-0.6732)*np.power(1.0+np.sqrt(T/3.583E+07),1.0+0.6732)) + np.power(T,-1.5)*(2.801E-03*np.exp(-4.198E+06/T)+4.362E-02*np.exp(-5.516E+06/T)+1.117E-03*np.exp(-8.050E+06/T))
		R87 = 6.170e-10/(np.sqrt(T/1.316E+02)*np.power(1.0+np.sqrt(T/1.316E+02),1.0-0.7481)*np.power(1.0+np.sqrt(T/3.427E+07),1.0+0.7481))

		tr = inputFile['trkev'].value
		tr = 1000.0*tr
		flux = 2.36e23*tr*tr*tr

		alpha21 = (NI  + NII) *R21  / (flux* 1.0 * 1.0e-18)
		alpha32 = (NII + NIII)*R32  / (flux* 0.8 * 1.0e-18)
		alpha43 = (NIII+ NIV) *R43  / (flux* 0.6 * 1.0e-18)
		alpha54 = (NIV + NV)  *R54  / (flux* 0.23 * 1.0e-18)
		alpha65 = (NV  + NVI) *R65  / (flux* 0.12 * 1.0e-18)
		alpha76 = (NVI + NVII)*R76  / (flux* 0.015 * 1.0e-18)
		alpha87 = (NVII+ NVIII)*R87 / (flux* 0.0027 * 1.0e-18)

		if field=='alpha21':
			return alpha21
		if field=='alpha32':
			return alpha32
		if field=='alpha43':
			return alpha43
		if field=='alpha54':
			return alpha54
		if field=='alpha65':
			return alpha65
		if field=='alpha76':
			return alpha76
		if field=='alpha87':
			return alpha87
	elif field=='gradTr':
		#Flux ~ Er/dx
		#Er  = np.array(inputFile['erad'].value)
		Tr  = np.array(inputFile['trkev'].value)*1000.
		dx  = np.array(inputFile['dx'].value)
		#The 4.2e-3 comes from Mayer 1965 (for N at 5ev and 6.0e-3 gm/cm^3)
		#Flux = 4.2e-3*(Er/dx)/(2.7*Tr)
		Flux = (4.2e-3*7.5657e-16*(Tr/8.617e-5)**3*(Tr/(8.617e-5*dx)))/(2.7*Tr)
		return Flux
	else:
		data = inputFile[field].value
	return data

def makePlotData(inputFile,field,logField=False,clabel="",outputFileName="Try.png"):
	inputFile = h5py.File(inputFile,"r")

	xpos = getData(inputFile,"posx")
	nb,nx,ny,nz = xpos.shape
	xmin = np.min(xpos)
	xmax = np.max(xpos)
	ypos = getData(inputFile,"posy")
	ymin = np.min(ypos)
	ymax = np.max(ypos)
	var = getData(inputFile,field)
	xpos = xpos.reshape( nb*nx*ny*nz )
	ypos = ypos.reshape( nb*nx*ny*nz )
	var  = var.reshape(  nb*nx*ny*nz )

	ratio = (ymax-ymin)/(xmax-xmin)

	NY=512
	NX=int(NY/ratio)

	grid_x = np.linspace(xmin,xmax,NX)
	grid_y = np.linspace(ymin,ymax,NY)
	grid_x,grid_y = np.meshgrid(grid_x,grid_y,)
	grid_v = griddata( (xpos,ypos), var, (grid_x,grid_y), fill_value=1.0e-10 )
	plotPS(grid_x,grid_y,grid_v,log=logField,clabel=clabel,outputFileName=outputFileName)

	shape = grid_x.shape
	shape=shape[0]*shape[1]
	grid_x = grid_x.reshape(shape)
	grid_y = grid_y.reshape(shape)
	grid_v = grid_v.reshape(shape)
	data = np.zeros((3,len(grid_x)))
	data[0,:]=grid_x
	data[1,:]=grid_y
	data[2,:]=grid_v
	ln = outputFileName[:-4]+'_2d.txt'
	np.savetxt( ln, data.T)

def makeLineData(inputFile,field,logField=False,clabel="",outputFileName="Try.png"):
	inputFile = h5py.File(inputFile,"r")
	xpos = getData(inputFile,"posx")
	nb,nx,ny,nz = xpos.shape
	xmin = np.min(xpos)
	xmax = np.max(xpos)
	ypos = getData(inputFile,"posy")
	ymin = np.min(ypos)
	ymax = np.max(ypos)
	var = getData(inputFile,field)
	xpos = xpos.reshape( nb*nx*ny*nz )
	ypos = ypos.reshape( nb*nx*ny*nz )
	var = var.reshape( nb*nx*ny*nz )

	ratio = (ymax-ymin)/(xmax-xmin)
	NY=512
	NX=int(NY/ratio)

	xmax = 1000.0

	grid_x = np.linspace(xmin,xmax,NX)
	grid_y = 250
	grid_v = griddata( (xpos,ypos), var, (grid_x,grid_y), fill_value=1.0e-10 )

# 	if field == 'optical_depth':
# 		ngrid = np.zeros(len(grid_v))
# 		dx    = np.zeros(len(grid_v))
# 		for i in xrange(len(ngrid)-1):
# 			dx[i] = grid_x[i+1]-grid_x[i]
# 		for i in xrange(len(ngrid)):
# 			ngrid[i] = np.sum(grid_v[0:i]*dx[0:i])
# 		ngrid = ngrid * 0.5 * 1.0e-18 * 1.0e-4
# 		grid_v = ngrid



	plotLine(grid_x,grid_y,grid_v,logField=logField,clabel=clabel,outputFileName=outputFileName)
	ln = outputFileName[:-4]+'.txt'
	data = np.zeros((2,len(grid_x)))
	data[0,:]=grid_x
	data[1,:]=grid_v
	np.savetxt( ln, data.T)

def plotLine(x,y,v,logField=False,clabel="Variable",outputFileName="Try.png"):
	if logField == True:
		v = np.log10(v)
	fig = plt.figure()
	plt.clf()
	plt.plot(x,v)
	plt.xlabel("X [um]")
	plt.ylabel(clabel)
	print "Saving: ", outputFileName
	plt.savefig(outputFileName)


def plotPS(x,y,v,log=False,clabel="Variable",outputFileName="Try.png"):
	if log == True:
		v = np.log10(v)
	nx,ny = np.shape(v)
	ratio=float(nx)/float(ny)
	fig = plt.figure()
	plt.clf()
	ax=plt.gca()
	ax.set_aspect("equal",adjustable="box")
	vmin = np.min(v)
	vmax = np.max(v)
	plt.pcolormesh(x,y,v,vmin=vmin,vmax=vmax)
	cbar = plt.colorbar(ax=ax,shrink=ratio)
	cbar.ax.set_ylabel(clabel)
	plt.xlabel("X [um]")
	plt.ylabel("Y [um]")
	print "Saving: ", outputFileName
	plt.savefig(outputFileName)



def main():
	batlFiles = glob.glob("z=0_var_1_n*.batl")
	batlFiles.sort()
	batlFiles = batlFiles[::5]

	#Vars = np.array(["flux","front_velocity","rho","trkev","zavg","level"])
	#Vars = np.array(["alpha"])
	#Vars = np.array(["gradTr","erad"])
	#Vars = np.array(["flux"])
	#Varl = np.array([True,False,False,False,False,False])
	#Varl = np.array([False,False])
	#Varl = np.array([False])

	Vars = np.array(["alpha21","alpha32","alpha43","alpha54","alpha65","alpha76","alpha87"])
	Varl = np.array([False,False,False,False,False,False,False])


	for i in xrange(len(Vars)):
		for j in xrange(len(batlFiles)):
			ln = "lineouts/"+Vars[i]+"_lines_%i.png"%(j)
			makeLineData(batlFiles[j],Vars[i],logField=Varl[i],clabel=Vars[i],outputFileName=ln)




# 	Vars = np.array(["rho","trkev","zavg"])
# 	Varl = np.array([False,False,False])



	# # Vars = np.array(["level"])
	# # Varl = np.array([False])

# 	for i in xrange(len(Vars)):
# 		for j in xrange(len(batlFiles)):
# 			ln = Vars[i]+"_lines_%i.png"%(j)
# 			makeLineData(batlFiles[j],Vars[i],logField=Varl[i],clabel=Varl[i],outputFileName=ln)

# 	batlFiles = "z=0_var_1_n00000000.batl"
# 	Vars = "rho"
# 	Varl = False
# 	makeLineData(batlFiles,Vars,logField=False,clabel="Rho",outputFileName="Try.png")

	#batlFiles = glob.glob("z=0_var_1_n*.batl")
	#batlFiles.sort()
# 	Vars = np.array(["rho","level","zavg","trkev"])
# 	VarL = np.array([True,False,False,False])
#
# 	for i in xrange(len(Vars)):
# 		for j in xrange(len(batlFiles)):
# 			ln = Vars[i]+"_%i.png"%(j)
# 			makePlotData(batlFiles[j],Vars[i],logField=VarL[i],clabel=Vars[i],outputFileName=ln)






if __name__== "__main__":
	main()