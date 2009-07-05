# -*- coding: utf-8 -*-
import pylab,numpy,math
from Scientific.IO import NetCDF


ncFile   = NetCDF.NetCDFFile('./simulation.nc', 'r')

f = pylab.figure(4,figsize=(22,4.5), dpi=100)

image = f.add_subplot(111)
propagation = ncFile.variables['Impulse']
propagationData = propagation.getValue()
	
nf = propagationData.shape[0]
nr = propagationData.shape[1]
nz = propagationData.shape[2]
	
data = numpy.zeros((nf/2,nz),dtype=float)

data = numpy.transpose(numpy.sum(propagationData,2))

print propagationData	


colormap = image.imshow(data,interpolation='bilinear',cmap=pylab.cm.jet,aspect='auto')
image.set_title('Should increase from left to right')

f.savefig('netcdftest.png')
