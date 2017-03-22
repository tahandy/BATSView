#!/usr/bin/env python3
from BATL_idl import BATL_idl
from BATL_hdf import BATL_hdf
import numpy as np
import matplotlib.pyplot as plt

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
	# cbar = plt.colorbar(ax=ax,shrink=ratio,orientation='horizontal')
	# cbar.ax.set_ylabel(clabel)
	plt.xlabel("X [um]")
	plt.ylabel("Y [um]")
	print("Saving: {:s}".format(outputFileName))
	plt.savefig(outputFileName,dpi=600)


# obj = BATL_hdf('../z=0_var_1_n00000000.batl',['rho','level'])


obj = BATL_idl('../2d_idl.out',['rho','level'])
X,Y,Z,D = obj.getUniform('level',nPts=[400,200,1])
plotPS(X,Y,D,clabel='level')
