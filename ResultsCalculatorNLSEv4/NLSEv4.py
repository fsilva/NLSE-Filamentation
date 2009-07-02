
import pytz
import matplotlib   
import pylab
import numpy
import scipy
import cairo
import sys
import os
import math
import commands
from Scientific.IO import NetCDF
import time






def Process(directory,name):
	try:
		f = open("../Results/"+directory+".png",'r')
		f.close()
		print 'Directory ' + directory + ' already processed.'
		return
	except:
		print 'Processing directory: ' + directory

	if(1):
#	try:
		try:
			simulationFile   = NetCDF.NetCDFFile(os.getcwd() + '/' + directory + '/simulation.nc', 'r')
		except:
			print 'Error in simulation files - sim probably not done yet'
			return

		try:
			finished    = getattr(simulationFile, 'finished') 
		except:
			print "Simulation not done yet"	
			return

#		pylab.clf()

		files = []
		drawImpulses(simulationFile,'impulses.png')
		files.append('impulses.png')
		drawSpectra(simulationFile,'spectra.png')
		files.append('spectra.png')
		#fwhmVal = [10,20]#drawFWHM(simulationFile,'fwhm.png')
		#files.append('fwhm.png')
		#drawStepAndError(simulationFile,'stepError.png')
		#files.append('stepError.png')
		drawPropagation(simulationFile,'impulse.png','spatialFWHM.png')
		files.append('impulse.png')
		files.append('spatialFWHM.png')
		drawCenterSpectrum(simulationFile,'cspectrum.png')
		files.append('cspectrum.png')
		drawCenterTimeProfile(simulationFile,'ctimeprofile.png')
		files.append('ctimeprofile.png')
		drawEnergy(simulationFile,'energy.png')
		files.append('energy.png')
		drawSpatialProfile(simulationFile,'profile.png')
		files.append('profile.png')
		drawPeakIntensity(simulationFile,'peakIntensity.png')
		files.append('peakIntensity.png')
		drawPeakElectronDensity(simulationFile,'peakRho1.png','peakRho2.png')
		files.append('peakRho1.png')
		files.append('peakRho2.png')


		phi = getPhi(simulationFile,0)

		version    = getattr(simulationFile, 'version') 
		maxPhi    = [01101011]#getattr(simulationFile, 'maxPhi') 
		deltaPhi  = [12312]#getattr(simulationFile, 'deltaPhi')
		timescale  = getattr(simulationFile, 'timescale')
		distancescale  = getattr(simulationFile, 'distancescale')
		deltaT    = getattr(simulationFile, 'deltaT') 
		zDistance  = getattr(simulationFile, 'zDistance')
		GVD    = getattr(simulationFile, 'GVD') 
		n2  = getattr(simulationFile, 'n2')
		absorption = getattr(simulationFile, 'absorption')
		lambda0  = getattr(simulationFile, 'lambda0')
		beamFWHM  = getattr(simulationFile, 'beamFWHM')
		boundaryRatio  = getattr(simulationFile, 'boundaryRatio')
		boundarySlope  = getattr(simulationFile, 'boundarySlope')
		boundaryI0     = getattr(simulationFile, 'boundaryI0')
		n               = getattr(simulationFile, 'n')
		sigma           = getattr(simulationFile, 'sigma')
		Ui              = getattr(simulationFile, 'Ui')
		sigmaK          = getattr(simulationFile, 'sigmaK')
		K               = getattr(simulationFile, 'K')
		rho_neutral     = getattr(simulationFile, 'rho_neutral')
		tau_r           = getattr(simulationFile, 'tau_r')
		tau_c           = getattr(simulationFile, 'tau_c')
		p 	        = getattr(simulationFile, 'p')


		print 'Rendering...'
		plot1 = cairo.ImageSurface.create_from_png("impulse.png")
		plot2 = cairo.ImageSurface.create_from_png("ctimeprofile.png")
		plot3 = cairo.ImageSurface.create_from_png("cspectrum.png")
		plot4 = cairo.ImageSurface.create_from_png("energy.png")
		plot5 = cairo.ImageSurface.create_from_png("profile.png")
		plot6 = cairo.ImageSurface.create_from_png("spatialFWHM.png")
		plot7 = cairo.ImageSurface.create_from_png("peakIntensity.png")
		plot8 = cairo.ImageSurface.create_from_png("impulses.png")
		plot9 = cairo.ImageSurface.create_from_png("spectra.png")
		plot10= cairo.ImageSurface.create_from_png("peakRho1.png")
		plot11= cairo.ImageSurface.create_from_png("peakRho2.png")

		#plot2 = cairo.ImageSurface.create_from_png("spectra.png")
		#plot3 = cairo.ImageSurface.create_from_png("fwhm.png")
		#plot4 = cairo.ImageSurface.create_from_png("stepError.png")
		#plot5 = cairo.ImageSurface.create_from_png("impulse.png")
		#plot6 = cairo.ImageSurface.create_from_png("spectrum.png")
		


		W = 4000.
		H = 4000*10./16. #FIXME

		GraphsX0 = 300.

		MapsY = H/2.


		surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, int(W), int(H))
		ctx = cairo.Context(surface)
		ctx.set_source_rgb(1, 1, 1)
	        ctx.rectangle(0, 0, W, H)
	        ctx.fill()

		x = GraphsX1 = GraphsX0
		y = 0
#plots
#Peak Rho1
		width  = plot10.get_width()
		GraphsY0 = height = plot10.get_height()
		ctx.save()
		ctx.set_source_surface(plot10,x,y)
		ctx.paint()
		ctx.restore()
		GraphsX1 += width*0.95
		y += height

#spatial profile
		width  = plot1.get_width()
		height = plot1.get_height()
		ctx.save()
		ctx.set_source_surface(plot1,x,y)
		ctx.paint()
		ctx.restore()
		y += height
	
#time profile	
		width  = plot2.get_width()
		height = plot2.get_height()
		ctx.save()
		ctx.set_source_surface(plot2,x,y)
		ctx.paint()
		ctx.restore()
		y += height
#spectral profile		
		height = plot3.get_height()
		ctx.save()
		ctx.set_source_surface(plot3,x,y)
		ctx.paint()
		ctx.restore()

# plots f(z) (on the right)

		y = 0
		x = GraphsX1

#power
		width  = plot4.get_width()
		height = plot4.get_height()
		ctx.save()
		ctx.set_source_surface(plot4,x,y)
		ctx.paint()
		ctx.restore()
		y += height

#Spatial FWHM
		width  = plot6.get_width()
		height = plot6.get_height()
		ctx.save()
		ctx.set_source_surface(plot6,x,y)
		ctx.paint()
		ctx.restore()
		y += height

#Peak Power
		width  = plot7.get_width()
		height = plot7.get_height()
		ctx.save()
		ctx.set_source_surface(plot7,x,y)
		ctx.paint()
		ctx.restore()
		y += height

#Peak Rho2
		width  = plot11.get_width()
		height = plot11.get_height()
		ctx.save()
		ctx.set_source_surface(plot11,x,y)
		ctx.paint()
		ctx.restore()
		y += height

#plots f(r), f(t) (bottom)
		x = GraphsX0
		y = GraphsY0*4

#temporal shape
		width  = plot8.get_width()
		height = plot8.get_height()
		ctx.save()
		ctx.set_source_surface(plot8,x,y)
		ctx.paint()
		ctx.restore()
		x += width

#spectra
		width  = plot9.get_width()
		height = plot9.get_height()
		ctx.save()
		ctx.set_source_surface(plot9,x,y)
		ctx.paint()
		ctx.restore()
		x += width
		
#Spatial Profile
		width  = plot5.get_width()
		height = plot5.get_height()
		ctx.save()
		ctx.set_source_surface(plot5,x,y)
		ctx.paint()
		ctx.restore()


#text
		ctx.scale(1,1)
		ctx.set_source_rgba(0, 0, 0, 1)
		ctx.select_font_face("FreeSans",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_BOLD)

		text = "Nonlinear Impulse Propagation Simulation - " + name
		extents = ctx.text_extents(text)
		
		ctx.set_font_size(20.0);
		ctx.move_to(W/2-extents[2], 20)
		ctx.show_text(text);

		ctx.select_font_face("FreeSans",cairo.FONT_SLANT_NORMAL,cairo.FONT_WEIGHT_NORMAL)

		ctx.set_font_size(24.0);


		strings1 = []
		strings2 = []
		strings1.append(version)
		strings2.append('')
		strings1.append('')
		strings2.append('')
		strings1.append("dT           ")
		strings2.append("= "+str(deltaT[0])+timescale)
		strings1.append("impulseWidth ")
		strings2.append("= "+str(deltaPhi[0])+timescale)
		strings1.append("zDistance    ")
		strings2.append("= "+str(zDistance[0])+distancescale)
		strings1.append("GVD          ")
		strings2.append("= "+str(GVD[0])+timescale+"^2"+distancescale+"^-1")
		strings1.append("n2           ")
		strings2.append("= "+str(n2[0]))
#		strings1.append("gamma         ")
#		strings2.append("= "+str(gamma[0]))
		#strings1.append("               "+distancescale+"^4"+timescale+"^-3KgA^-1")
		strings1.append("absorption   ")
		strings2.append("= "+str(absorption[0])+distancescale+'^-1')
		strings1.append("lambda0       ")
		strings2.append("= "+str(lambda0[0])+timescale+'^-1')
		strings1.append("beamFWHM    ")
		strings2.append("= "+str(beamFWHM[0])+distancescale)
		strings1.append('')
		strings2.append('')
		strings1.append("boundaryRatio    ")
		strings2.append("= "+str(boundaryRatio[0]))
		strings1.append("boundarySlope    ")
		strings2.append("= "+str(boundarySlope[0]))
		strings1.append("boundaryI0    ")
		strings2.append("= "+str(boundaryI0[0]))
		strings1.append('')
		strings2.append('')
		#strings1.append("FWHM (z=0)")
		#strings2.append("= %3.3e"%fwhmVal[0]+timescale)
		#strings1.append("FWHM (z=Max)")
		#strings2.append("= %3.3e"%fwhmVal[1]+timescale)
		strings1.append('')
		strings2.append('')
		#strings1.append("phi (z=L,t=0)")
		#strings2.append("= %3.3fpi"%(phi/numpy.pi))


		strings1.append('n')
		strings2.append("= "+str(n[0]))
		strings1.append('sigma')
		strings2.append("= "+str(sigma[0])+" " + distancescale + "^2")
		strings1.append('Ui')
		strings2.append("= "+str(Ui[0])+" J")
		strings1.append('sigmaK')
		strings2.append("= "+str(sigmaK[0])+" " + distancescale + "^2")
		strings1.append('K')
		strings2.append("= "+str(K[0]))
		strings1.append('rho_neutral')
		strings2.append("= "+str(rho_neutral[0])+" " + distancescale + "^-3")
		strings1.append('tau_r')
		strings2.append("= "+str(tau_r[0])+" " + timescale)
		strings1.append('tau_c')
		strings2.append("= "+str(tau_c[0])+" " + timescale)
		strings1.append('p')
		strings2.append("= "+str(p[0]))


		for i in range(len(strings1)):
			ctx.move_to(5, 90+i*30)
			ctx.show_text(strings1[i]);

		for i in range(len(strings2)):
			ctx.move_to(185, 90+i*30)
			ctx.show_text(strings2[i]);
	
		max_y = 90+len(strings2)*20



#		surface.write_to_png("./Results/"+directory+".png")
		#surface.write_to_png("./Results/"+directory+".png")
		surface.write_to_png("../Results/"+directory+".png")
		surface.write_to_png("./Results/"+directory+".png")

		for f in files:
			os.system('rm '+f)
#	except:
#		print "Error processing directory ", directory


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################


def drawImpulses(ncFile,outfile):
	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	#print propagationData.shape

	nt = propagationData.shape[0]

	profile_i = numpy.zeros((nt/2),dtype=float)
	profile_f = numpy.zeros((nt/2),dtype=float)

	c = 3e8
	epsilon0 = 8.85e-12

	for k in xrange(nt/2):
		real = propagationData[k*2+0][0][0]
		imag = propagationData[k*2+1][0][0]
		profile_i[k] = math.sqrt(real**2+imag**2)
		real = propagationData[k*2+0][0][-1]
		imag = propagationData[k*2+1][0][-1]
		profile_f[k] = math.sqrt(real**2+imag**2)

	f = pylab.figure(1,figsize=(6,4.5),dpi=100)
	f.clf()

	plot = f.add_subplot(111)   #this nt is 2x the simulation nt (complex numbers)
	plot.plot(numpy.arange(-nt/4,+nt/4,1)*float(deltaT),profile_i,c='b')
	plot.plot(numpy.arange(-nt/4,+nt/4,1)*float(deltaT),profile_f,c='r')
	plot.set_ylabel('E (V' + distancescale+'^-1)')
	plot.set_xlabel('Time ('+timescale+')')
	plot.set_title('Initial and Final Temporal Shapes at r=0')

	f.savefig(outfile)
	del plot,f


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def drawSpectra(ncFile,outfile):
	propagation = ncFile.variables['Spectrum']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	#print propagationData.shape

	nf = propagationData.shape[0]

	profile_i = numpy.zeros((nf/2),dtype=float)
	profile_f = numpy.zeros((nf/2),dtype=float)

	c = 3e8
	epsilon0 = 8.85e-12

	for k in xrange(nf/2):
		real = propagationData[k*2+0][0][0]
		imag = propagationData[k*2+1][0][0]
		if(k >= nf/4):
			profile_i[k-nf/4] = math.sqrt(real**2+imag**2)
		else:
			profile_i[k+nf/4] = math.sqrt(real**2+imag**2)
		real = propagationData[k*2+0][0][-1]
		imag = propagationData[k*2+1][0][-1]
		if(k >= nf/4):
			profile_f[k-nf/4] = math.sqrt(real**2+imag**2)
		else:
			profile_f[k+nf/4] = math.sqrt(real**2+imag**2)

	f = pylab.figure(1,figsize=(6,4.5),dpi=100)
	f.clf()

	plot = f.add_subplot(111)   #this nt is 2x the simulation nt (complex numbers)
	plot.plot(numpy.arange(-nf/4,+nf/4,1)/float(deltaT),profile_i,c='b')
	plot.plot(numpy.arange(-nf/4,+nf/4,1)/float(deltaT),profile_f,c='r')
	plot.set_ylabel('E (V' + distancescale+'^-1)')
	plot.set_xlabel('Time ('+timescale+'^-1)')
	plot.set_title('Initial and Final Spectra at r=0')

	f.savefig(outfile)
	del plot,f


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################


def drawFWHM(ncFile, outfile):
	f = pylab.figure(2,figsize=(6,4.5), dpi=100)
	f.clf()

	deltaT         = getattr(ncFile, 'deltaT') 
	timescale      = getattr(ncFile, 'timescale')
	plot = f.add_subplot(111)
	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	dataR = propagationData[0::2]
	dataI = propagationData[1::2]
	nx = dataR.shape[0]
	max_ = numpy.zeros((dataR.shape[1]))
	min_ = numpy.zeros((dataR.shape[1]))
	fwhm = numpy.zeros((dataR.shape[1]))
	fwhm2 = numpy.zeros((dataR.shape[1]))
	for j in xrange(dataR.shape[1]):
		max__ = dataR[0][j]
		min__ = dataR[0][j]
		for i in xrange(dataR.shape[0]):
			dataR[i][j] = numpy.sqrt(dataR[i][j]**2+dataI[i][j]**2)
			if(dataR[i][j] > max__):
				max__ = dataR[i][j]
			if(dataR[i][j] < min__):
				min__ = dataR[i][j]
		min_[j] = min__
		max_[j] = max__

	for j in xrange(dataR.shape[1]):
		fwhm[j] = 0
		fwhm2[j] = 0
		i0 = -1
		i1 = -1
		for i in xrange(dataR.shape[0]):
			if(dataR[i][j]-min_[j] >= 0.5*(max_[j]-min_[j])):
				fwhm2[j] += 1
			if(dataR[i][j]-min_[j] >= 0.5*(max_[j]-min_[j]) and i0 == -1):
				i0 = i
			if(dataR[i][j]-min_[j] <= 0.5*(max_[j]-min_[j]) and i1 == -1 and i0 != -1):
				i1 = i
		j0 = i0+(0.5*(max_[j]-min_[j])-dataR[i0][j])/(dataR[i0][j]-dataR[i0-1][j])
		j1 = i1+(0.5*(max_[j]-min_[j])-dataR[i1][j])/(dataR[i1][j]-dataR[i1-1][j])
			
		fwhm[j] = j1-j0

	deltaT = getattr(ncFile, 'deltaT') 

#	plot.clf()
	plot.plot(numpy.arange(len(fwhm))*deltaT,fwhm*deltaT[0],c='b')
	plot.set_ylabel('FWHM ('+timescale+')')
	plot.set_xlabel('Time ('+timescale+')')
	plot.set_title('FWHM vs Time')
#	pylab.plot(fwhm2*deltaT[0],c='r')
#	print 'FWHM - initial = ',fwhm[0]*deltaT[0]*0.001 ," final = ",fwhm[-1]*deltaT[0]*0.001

	f.savefig(outfile)
	del plot
	del f
	return [fwhm[0]*deltaT[0],fwhm[-1]*deltaT[0]]


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawStepAndError(ncFile,outfile):
	f = pylab.figure(3,figsize=(6,9), dpi=100)
	f.clf()

	deltaT         = getattr(ncFile, 'deltaT') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')

	plot1 = f.add_subplot(211)
	plot1.set_title('Step and Error vs Time')

	propagation = ncFile.variables['zStep']
	data = propagation.getValue()
	plot1.plot(numpy.arange(len(data))*deltaT,data)
	plot1.set_ylabel('Step ('+distancescale+')')
	plot1.set_xlabel('Time ('+timescale+')')

	plot2 = f.add_subplot(212)
	propagation = ncFile.variables['Error']
	data = propagation.getValue()
	plot2.plot(numpy.arange(len(data))*deltaT,data)
	plot2.set_ylabel('Error')
	plot2.set_xlabel('Time ('+timescale+')')
#	plot2.set_title('Power vs Time')

	f.savefig(outfile)

	del plot1,plot2,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawCenterSpectrum(ncFile,outfile):
	f = pylab.figure(4,figsize=(22,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	propagation = ncFile.variables['Spectrum']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')
	
	nf = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]
	
	data = numpy.zeros((nf/2,nz),dtype=float)
	
	for i in xrange(nz):
		for k in xrange(nf/2):
			real = propagationData[k*2+0][0][i]
			imag = propagationData[k*2+1][0][i]
			if(k >= nf/4):
				data[k-nf/4][i] = math.sqrt(real*real+imag*imag)#numpy.log(math.sqrt(real*real+imag*imag))
			else:
				data[nf/4+k][i] = math.sqrt(real*real+imag*imag)#numpy.log(math.sqrt(real*real+imag*imag))
	
	

	#dataR = propagationData[0::2]
	#dataI = propagationData[1::2]
	#data  = numpy.zeros(dataR.shape)
	#nx = dataR.shape[0]
	#for i in xrange(data.shape[0]):
		#for j in xrange(data.shape[1]):
			#if(i >= nx/2):
				#data[i-nx/2][j] = numpy.sqrt(dataR[i][j]**2+dataI[i][j]**2)
			#else:
				#data[nx/2+i][j] = numpy.sqrt(dataR[i][j]**2+dataI[i][j]**2)
	#absVal = scipy.misc.pilutil.imresize(data,(nx,nx*3))
	freq = 2.*numpy.pi/2./float(deltaT)
	colormap = image.imshow(data,interpolation='bilinear',cmap=pylab.cm.jet,extent=(0, float(distance), -freq,+freq),aspect='auto')
	image.set_ylabel(r"$\omega$ ("+timescale+'^-1)')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Impulse Propagation (frequency domain)')

	#f.colorbar(colormap)
	
	f.savefig(outfile)

	del image,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawPropagation(ncFile,outfile,outfile2):
	f = pylab.figure(5,figsize=(22,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	#print propagationData.shape

	nt = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]

	#print nt,nr,nz

	data = numpy.zeros((nr*2,nz),dtype=float)

#	print numpy.sum(propagationData[0::2][0][0]**2+propagationData[1::2][0][0]**2)
	#print propagationData[0].shape,propagationData[0][0][0].shape

	c = 3e8
	epsilon0 = 8.85e-12

	for i in xrange(nz):
		for j in xrange(nr):
			energy = 0
			for k in xrange(nt/2):
				real = propagationData[k*2+0][j][i]
				imag = propagationData[k*2+1][j][i]
				energy += c*epsilon0*(real*real+imag*imag)*deltaT[0]*deltaR[0]
			#print energy," ",
			data[nr-j-1][i] = data[nr+j][i] = energy#numpy.log(energy)
	#dataR = propagationData[0::2]
	#dataI = propagationData[1::2]
	#nx = dataR.shape[0]
	#for i in xrange(dataR.shape[0]):
	#	for j in xrange(dataR.shape[1]):
	#		dataR[i][j] = numpy.sqrt(dataR[i][j]**2+dataI[i][j]**2)
#	absVal = scipy.misc.pilutil.imresize(dataR,(nx,nx*3))


	#time = len(dataR)/2.*float(deltaT)
	colormap = image.imshow(data,interpolation='bilinear',cmap=pylab.cm.jet,extent=(0, float(distance[0]), -nr*deltaR[0],+nr*deltaR[0]),aspect='auto')
	image.set_ylabel('r ('+distancescale+')')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Impulse Propagation (spacial domain)')

	#f.colorbar(colormap)
	
	f.savefig(outfile)

	##################33
	#### calc FWHM

	max_ = numpy.zeros((nz))
	min_ = numpy.zeros((nz))
	fwhm = numpy.zeros((nz))
	z = numpy.arange(0,distance[0],distance[0]/nz)

	#ff = open('blah','w')
	#for i in xrange(nr):
#		ff.write(str(data[i][0])+'\n')
#	ff.close()

	for j in xrange(nz):
		max__ = data[0][j]
		min__ = data[0][j]
		for i in xrange(nr):
			if(data[i][j] > max__):
				max__ = data[i][j]
			if(data[i][j] < min__):
				min__ = data[i][j]
		min_[j] = min__
		max_[j] = max__
	
	for j in xrange(nz):
		fwhm[j] = 0
		i0 = -1
		for i in xrange(nr):
			if(data[nr+i][j]-min_[j] <= 0.5*(max_[j]-min_[j]) and i0 == -1):
				i0 = i
		j0 = i0+(0.5*(max_[j]-min_[j])-data[nr+i0][j])/(data[nr+i0][j]-data[nr+i0-1][j])

		fwhm[j] = 2*(j0+1)*deltaR[0]

		print i0,j0,'        ',data[nr+i0][j]/max_[j],data[nr+i0-1][j]/max_[j], '   z=',distance[0]/nz*j,'  fwhm=',fwhm[j]

	f = pylab.figure(9,figsize=(11,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	image.plot(z,fwhm)
	image.set_ylabel('FWHM ('+distancescale+')')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Impulse FWHM (spacial domain)')

	f.savefig(outfile2)
	del image,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawPeakIntensity(ncFile,outfile):

	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	nt = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]

	peakIntensity = numpy.zeros((nz))
	z = numpy.arange(0,distance[0],distance[0]/nz)

	c = 3e8
	print 'TODO: substitiuir pelo n do argon'
	n = 1
	epslon0 = 8.85418781e-12

	for i in xrange(nz):
		peakIntensity[i] = -1e20
		for j in xrange(nr):
			for k in xrange(nt/2):
				real = propagationData[k*2+0][j][i]
				imag = propagationData[k*2+1][j][i]
				#area = 2.0*3.1415*float(j+1)*deltaR[0]*deltaR[0]
				#I = c*n*epslon0/2*(real**2+imag**2)/area
				I = c*n*epslon0/2*(real**2+imag**2)
				if(I > peakIntensity[i]):
					peakIntensity[i] = I

	f = pylab.figure(11,figsize=(11,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	image.plot(z,peakIntensity)
	image.set_ylabel('Peak Intensity (W/m^2)')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Peak Temporal Intensity over distance')

	f.savefig(outfile)

	del image,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawSpatialProfile(ncFile,outfile):
	f = pylab.figure(8,figsize=(6,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	#print propagationData.shape

	nt = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]

	#print nt,nr,nz

	data_i = numpy.zeros((nr*2),dtype=float)
	data_f = numpy.zeros((nr*2),dtype=float)

#	print numpy.sum(propagationData[0::2][0][0]**2+propagationData[1::2][0][0]**2)
	#print propagationData[0].shape,propagationData[0][0][0].shape
	c = 3e8
	print 'TODO: substitiuir pelo n do argon'
	n = 1
	epslon0 = 8.85418781e-12

	for j in xrange(nr):
		amplitude = 0
		for k in xrange(nt/2):
			real = propagationData[k*2+0][j][-1]
			imag = propagationData[k*2+1][j][-1]
			I = c*n*epslon0/2*(real**2+imag**2)
			amplitude += I

		data_f[nr-j-1] = data_f[nr+j] = c*epslon0*n*amplitude/2
		amplitude = 0
		for k in xrange(nt/2):
			real = propagationData[k*2+0][j][0]
			imag = propagationData[k*2+1][j][0]
			I = c*n*epslon0/2*(real**2+imag**2)
			amplitude += I

		data_i[nr-j-1] = data_i[nr+j] = c*epslon0*n*amplitude/2

	image.plot(numpy.arange(-nr,nr,1)*deltaR,data_f,'r')
	image.plot(numpy.arange(-nr,nr,1)*deltaR,data_i,'b')
#	image.set_xlim((-nr/2*deltaR,+nr/2*deltaR))
	image.set_ylabel(r'$\sum_tIntensity(W/m^2)$')
	image.set_xlabel('r ('+distancescale+')')
	image.set_title('Spatial Profile at the end')
	
	f.savefig(outfile)

	del image,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################

def 	drawCenterTimeProfile(ncFile,outfile):
	f = pylab.figure(5,figsize=(22,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	propagation = ncFile.variables['Impulse']
	propagationData = propagation.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	#print propagationData.shape

	nt = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]

	#print nt,nr,nz

	data = numpy.zeros((nt/2,nz),dtype=float)

#	print numpy.sum(propagationData[0::2][0][0]**2+propagationData[1::2][0][0]**2)
	#print propagationData[0].shape,propagationData[0][0][0].shape

	for i in xrange(nz):
		for k in xrange(nt/2):
			real = propagationData[k*2+0][0][i]
			imag = propagationData[k*2+1][0][i]
			data[k][i] = math.sqrt(real*real+imag*imag)#numpy.log(math.sqrt(real*real+imag*imag))
			#print energy," ",
	#dataR = propagationData[0::2]
	#dataI = propagationData[1::2]
	#nx = dataR.shape[0]
	#for i in xrange(dataR.shape[0]):
	#	for j in xrange(dataR.shape[1]):
	#		dataR[i][j] = numpy.sqrt(dataR[i][j]**2+dataI[i][j]**2)
#	absVal = scipy.misc.pilutil.imresize(dataR,(nx,nx*3))

	#time = len(dataR)/2.*float(deltaT)
	colormap = image.imshow(data,interpolation='bilinear',cmap=pylab.cm.jet,extent=(0, float(distance), -nt/2*deltaT[0],+nt/2*deltaT[0]),aspect='auto')
	image.set_ylabel('t ('+timescale+')')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Impulse Propagation (r=0 time profile)')

	#f.colorbar(colormap)
	
	f.savefig(outfile)

	del image,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################


def drawEnergy(ncFile,outfile):
	propagation = ncFile.variables['Impulse']
	deltaT         = getattr(ncFile, 'deltaT') 
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	propagationData = propagation.getValue()

	nt = propagationData.shape[0]
	nr = propagationData.shape[1]
	nz = propagationData.shape[2]

	data = numpy.zeros((nz))
	for i in xrange(nz):
		data[i] = 0
		for k in xrange(nr):
			for j in xrange(nt/2):	
				real = propagationData[j*2+0][k][i]
				imag = propagationData[j*2+1][k][i]
				data[i] += k*(real*real+imag*imag)
		#print data[i]

	f = pylab.figure(8,figsize=(11,4.5),dpi=100)
	f.clf()
	plot = f.add_subplot(111)
	plot.plot(numpy.arange(len(data))*float(distance)/float(len(data)),data)
	plot.set_ylabel('Energy(J)')
	plot.set_xlabel('Distance('+distancescale+')')
	plot.set_title('Energy vs Distance')
	f.savefig(outfile)
	del plot,f

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################

def getPhi(ncFile,t):
	propagation = ncFile.variables['Impulse']
	deltaT         = getattr(ncFile, 'deltaT') 

	propagationData = propagation.getValue()
	dataR = propagationData[0::2]
	dataI = propagationData[1::2]


	maxtime = len(dataR)/2.*float(deltaT)
	if(t > maxtime or -t > maxtime):
		print "Error: phi asked at t outside the window - t=",t,"  maxT=",maxtime
		return 0
	i = int(t/maxtime)+len(dataR)/2

	#phi = numpy.angle(numpy.complex(dataR[i][-1],dataI[i][-1]))
	phi = 9999
	return phi
	

#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



def 	drawPeakElectronDensity(ncFile,outfile1,outfile2):

	rho = ncFile.variables['ElectronDensity']
	rhoData = rho.getValue()
	deltaT         = getattr(ncFile, 'deltaT') 
	deltaR         = getattr(ncFile, 'deltaR') 
	timescale      = getattr(ncFile, 'timescale')
	distancescale  = getattr(ncFile, 'distancescale')
	distance       = getattr(ncFile, 'zDistance')

	nt = rhoData.shape[0]
	nr = rhoData.shape[1]
	nz = rhoData.shape[2]

	peakRho1 = numpy.zeros((nz))
	z = numpy.arange(0,distance[0],distance[0]/nz)

	for i in xrange(nz):
		peakRho1[i] = 0
		for k in xrange(nt/2):
			P = numpy.log(rhoData[k][0][i])
			if(P > peakRho1[i]):
				peakRho1[i] = P

	f = pylab.figure(13,figsize=(11,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	image.plot(z,peakRho1)
	image.set_ylabel('Peak Electron Density (m^-3)')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Peak Electron Density over distance (@r=0)')

	f.savefig(outfile2)

# filament structure
	data = numpy.zeros((2*nr,nz),dtype=float)

	for i in xrange(nz):
		for j in xrange(nr):
			peak = 0
			for k in xrange(nt/2):
				if(peak < rhoData[k][j][i]):
					peak = rhoData[k][j][i]
			data[nr-j-1][i] = data[nr+j][i] = peak

	f = pylab.figure(12,figsize=(22,4.5), dpi=100)
	f.clf()

	image = f.add_subplot(111)
	image.imshow(data,interpolation='bilinear',cmap=pylab.cm.hot,extent=(0, float(distance), -nr*deltaR[0],+nr*deltaR[0]),aspect='auto')
	image.set_ylabel('r ('+distancescale+')')
	image.set_xlabel('Distance ('+distancescale+')')
	image.set_title('Peak Electron Density over distance')


	f.savefig(outfile1)

	del image,f


