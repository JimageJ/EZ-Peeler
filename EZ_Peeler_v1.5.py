"""
*******************************************************************************************************************************
											Epidermal Z stack Peeler (EZ Peeler)
															v 1.41
				
							Written by Jim Rowe whilst working in Stuart Casson's and Alexander Jones' labs
									(Uni of Sheffield and SLCU, Cambridge respectively)
											
						In case it's not obvious, I'm a biologist, not a coder, so apologies for inefficiencies 
													and sloppy spaghetti code!
										
														@BotanicalJim
													james.rowe@slcu.cam.ac.uk
										
														For references:
	Zoulias N, Brown J, Rowe J, Casson SA (2020) 'HY5 is not integral to light mediated stomatal development in Arabidopsis.'
	
					
					Development started 15/05/2017 when I should have been doing real physiology work.... 

********************************************************************************************************************************
"""

#import libraries required for plugin

from ij import IJ, ImagePlus, ImageStack, CompositeImage
from java.awt import Color
from ij.gui import Roi, PolygonRoi, NonBlockingGenericDialog, Overlay, ImageRoi
from ij.process import ImageProcessor, StackStatistics, ImageConverter, FloatProcessor,ColorProcessor, ByteProcessor
from ij.plugin import Slicer, ImageCalculator, Duplicator, ZProjector
from ij.plugin.filter import GaussianBlur, RankFilters, ThresholdToSelection
from array import array, zeros
from collections import OrderedDict
from script.imglib.math import Compute, Subtract, Divide, Multiply
from script.imglib import ImgLib  
from java.lang import Thread

from fiji.util.gui import GenericDialogPlus

try: 
	from net.haesleinhuepf.clij2 import CLIJ2
	from net.haesleinhuepf.clij import CLIJ
except:
	errorDialog("""EZ-Peeler requires clij to function. 
	
	To install please follow these instructions: 
	
	1. Click Help>Update> Manage update sites
	2. Make sure the "clij" and "clij2" update sites are selected.
	3. Click Close> Apply changes.
	4. Close and reopen ImageJ""")
clij2 = CLIJ2.getInstance()
clij = CLIJ.getInstance()
clij2.clear()

"""**********************************define functions*******************************************"""
def errorDialog(message):
	gd = NonBlockingGenericDialog("EZ Peeler - Error")

	gd.addMessage(message)
	gd.showDialog()
	return


def areamap(imp1, xscale, yscale):
	
	"""Calculates a areamap from the heightmap"""

	imp2= imp1.duplicate()
	imp3= imp1.duplicate()
	imp4= imp1.duplicate()
	imp5= imp1.duplicate()

	#work out change in height in several directions
	IJ.run(imp2, "Convolve...", "text1=[0\n1\n-1] ")
	IJ.run(imp3, "Convolve...", "text1=[0 0 0\n0 1 -1\n0 0 0] ")
	IJ.run(imp4, "Convolve...", "text1=[-1\n1\n0] ")
	IJ.run(imp5, "Convolve...", "text1=[0 0 0\n-1 1 0\n0 0 0] ")

	
	#push to graphics card
	simp2=clij2.push(imp2)
	simp3=clij2.push(imp3)
	simp4=clij2.push(imp4)
	simp5=clij2.push(imp5)
	
	dimp2=clij2.create(simp2)
	dimp3=clij2.create(simp3)
	dimp4=clij2.create(simp4)
	dimp5=clij2.create(simp5)

	# height^2 + scale^2 = edgelength^2
	
	# height^2
	clij2.power(simp2,dimp2, 2)
	clij2.power(simp3,dimp3, 2)
	clij2.power(simp4,dimp4, 2)
	clij2.power(simp5,dimp5, 2)

	# height^2 + scale^2
	clij2.addImageAndScalar(dimp2, simp2, yscale*yscale)
	clij2.addImageAndScalar(dimp3, simp3, xscale*xscale)
	clij2.addImageAndScalar(dimp4, simp4, yscale*yscale)
	clij2.addImageAndScalar(dimp5, simp5, xscale*xscale)

	# Calculate edgelength
	clij2.power(simp2,dimp2, 0.5)
	clij2.power(simp3,dimp3, 0.5)
	clij2.power(simp4,dimp4, 0.5)
	clij2.power(simp5,dimp5, 0.5)
	
	dimp6=clij2.create(dimp2)
	dimp7=clij2.create(dimp4)
	
	#trianglular polygon area = Xedgelength x Yedgelength/2 
	#
	clij2.multiplyImages(dimp2, dimp3, dimp6)
	clij2.multiplyImages(dimp4, dimp5, dimp7)
	clij2.addImages(dimp7, dimp6, dimp2)
	clij2.multiplyImageAndScalar(dimp2, dimp3, 0.5)

	#pull result from GPU
	areamapIMP=clij2.pull(dimp3)
	
	#cleanup, free GPU memory
	simp2.close()
	simp3.close()
	simp4.close()
	simp5.close()
	
	dimp2.close()
	dimp3.close()
	dimp4.close()
	dimp5.close()
	dimp6.close()
	dimp7.close()
	return areamapIMP



def extractChannel(imp, nChannel, nFrame, changeBitType):
	"""extract a channel from the image, returning a new 16 bit imagePlus labelled with the channel name"""

	stack = imp.getImageStack()
	ch=ImageStack(imp.width, imp.height)
	for i in range(imp.getNSlices()):
		index = imp.getStackIndex(nChannel, i, nFrame)
		ch.addSlice(str(i), stack.getProcessor(index))
	imp3 = ImagePlus("Channel " + str(nChannel), ch).duplicate()
	stats =StackStatistics(imp3) 
	IJ.setMinAndMax(imp3, stats.min, stats.max)
	if andOp == 0xFF:
		ImageConverter(imp3).convertToGray8()
	else:
		ImageConverter(imp3).convertToGray16() 	
	return imp3

def extractFrame(imp, nFrame):
	"""extract a frame from the image, returning a new 16 bit imagePlus labelled with the channel name"""

	stack = imp.getImageStack()
	fr=ImageStack(imp.width, imp.height)
	for i in range(1, imp.getNSlices() + 1):
		for nChannel in range(1, imp.getNChannels()+1):
			index = imp.getStackIndex(nChannel, i, nFrame)
			fr.addSlice(str(i), stack.getProcessor(index))
	imp3 = ImagePlus("Frame " + str(nFrame), fr).duplicate()
	imp3.setDimensions(imp.getNChannels(), imp.getNSlices(), 1)
	comp = CompositeImage(imp3, CompositeImage.COMPOSITE)  
	comp.show()
	return comp



def smoothFilter(imp, filtertype, radius):
	"""iterate through slices to perform either gaussian or median filter"""
	try:
		src = clij2.push(imp)
		dst = clij2.create(src)
	except:	
		try:
		
			Thread.sleep(500)
			print("Succeeded to sending to graphics card on the second time...")
			src = clij2.push(imp)
			dst = clij2.create(src)
		except:
			errorDialog("""Could not send image to graphics card, it may be too large!
		
			Easy solutions: Try	processing as 8-bit, cropping or scaling the image, or
			select a different CLIJ2 GPU.

			This issue is often intermittent, so trying again may also work! 

			See the "Big Images on x graphics cards' notes at:
			https://clij2.github.io/clij2-docs/troubleshooting for more solutions
			
			"""	+ str(clij2.reportMemory()) )
			"""

			put non clij filters here
			if filtertype == "Gaussian":
				clij2.blur(src, dst, radius, radius)
				#imp = clij2.pull(dst)
			if filtertype == "Median":
				clij2.medianBox(src, dst, radius, radius,0)
				#imp = clij2.pull(dst)
			if filtertype == "3D Median":
				clij2.medianBox(src, dst, radius, radius, radius)
				#imp = clij2.pull(dst)
			"""

	if filtertype == "2D Gaussian":
		clij2.blur2D(src, dst, radius, radius)
		#imp = clij2.pull(dst)
	if filtertype == "3D Gaussian":
		clij2.blur3D(src, dst, radius, radius, radius)
		#imp = clij2.pull(dst)
	if filtertype == "Median":
		clij2.median2DBox(src, dst, radius, radius)
		#imp = clij2.pull(dst)
	if filtertype == "3D Median":
		clij2.median3DBox(src, dst, radius, radius, radius)
		#imp = clij2.pull(dst)
	
	src.close()
	
	#return as imagePlus
	return dst
def findEdge(imp3, sobeltype):


	if sobeltype =="2D Sobel":
		IJ.run(imp3, "Find Edges", "")
	if sobeltype =="1D Sobel":
		IJ.run(imp3, "Convolve...", "text1=[1 2 1 \n0 0 0 \n-1 -2 -1\n] normalize stack")
	if sobeltype =="1 X 3 Gradient":
		IJ.run(imp3, "Convolve...", "text1=[1\n0\n-1] normalize stack")
	if sobeltype =="Laplace filter":
		IJ.run(imp3, "Convolve...", "text1=[-1 -1 -1\n-1 8 -1\n-1 -1 -1] normalize stack")
	if sobeltype =="1 X 7 Gradient":
		IJ.run(imp3, "Convolve...", "text1=[1\n1\n1\n0\n-1\n-1\n-1] normalize stack")
	if sobeltype =="3 X 7 Gradient":
		IJ.run(imp3, "Convolve...", "text1=[1 1 1\n1 1 1\n1 1 1\n0 0 0\n-1 -1 -1\n-1 -1 -1\n-1 -1 -1] normalize stack")
	if sobeltype =="1D mexican hat":
		IJ.run(imp3, "Convolve...", "text1=[-1 \n -2 \n 1 \n 4 \n 1 \n -2 \n -1] normalize stack")
	return imp3
	
def getOptions1(imp):
	"""Get user defined options for image preprocessing"""
                      
	gd = NonBlockingGenericDialog("EZ Peeler - Preprocessing")

	gd.addMessage("""                EZ Peeler is written by Jim Rowe
	
					These options define which channel and timeframe to use for segmentation
					and any preproceessing steps.
					
					Channel number should be the cell wall stain/membrane marker that best
					defines the surface.

					The Timepoint selects which frame to segment from timeseries data. After
					segmentating this frame, you will then have the option to apply these
					settings to segment the whole time series.

					Processing as an 8 bit stack allows faster segmentation and and overcomes
					memory limitations when dealing with large image sets (Recommended on).
					If switched off, stacks will be processed as 16-bit images.
					
					Filters are used to improve segmentation, and are not applied to the 
					orignal image. Smoothing filters work well at removing granular noise
					from confocal stacks. 

					Edge detection filters make surface detection easier in stacks of 
					non-uniform intensity but grandular noise can be made
					worse with these filters, so use with a median filter is 
					advisible.""")

	
	#create a list of the channels in the provided imagePlus
	types = []
	for i in xrange(imp.getNChannels()):
		types.append(str(i))
	#user can pick which channel to base the peeling on
		
	gd.addChoice("Channel number to use for segmentation", types, types[-1])
	#user can pick which timepoint to segment
	
	timepoint=[]
	for i in xrange(imp.getNFrames()):
		timepoint.append(str(i))

	gd.addChoice("Timepoint to segment", timepoint, timepoint[0])

	
	#user specifies the smoothing process and the size of the smoothing
	gd.addCheckbox("Process as 8-bit image", True)
	smoothtype = ["2D Gaussian", "3D Gaussian", "3D Median", "none"]
	gd.addChoice("Image smoothing type", smoothtype, smoothtype[0])
	
	smoothing=["1", "2", "3", "4", "5", "6", "7", "8", "9"]
	gd.addChoice("Image smoothing sigma", smoothing, smoothing[3])
	#Use specifies whether they want to use edge detection
	sobel=["none", "1D Sobel", "2D Sobel", "Laplace filter","1 X 3 Gradient", "1 X 7 Gradient", "3 X 7 Gradient", "1D mexican hat"]
	gd.addChoice("Edge detection filter", sobel, sobel[0])
	GPUs=clij.getAvailableDeviceNames()
	gd.addChoice("CLIJ2 GPU choice", GPUs, GPUs[0])
	 
	gd.showDialog()
	

	channel = gd.getNextChoice()
	changeBitType= gd.getNextBoolean()
	frame= gd.getNextChoice()
	filtertype = gd.getNextChoice()
	smoothsize = gd.getNextChoice()
	sobeltype = gd.getNextChoice()
	GPU = gd.getNextChoice()
	#assign bitwise operator to convert signed to unsigned stack image values.
	if changeBitType==1:
		andOp = 0xFF
	else:
		andOp = 0xffff
	# if dialog canceled then return nothing
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 4"
	else:
		canceled = 0
	if gd.wasOKed():
		oked=1
	else:
		oked=0
	return channel, frame, filtertype, smoothsize, canceled, sobeltype, andOp, GPU, oked


def getOptions2(imp, maxPixel,height, stack,andOp, width,defaultThreshold):


	"""Get user defined options for segmentation"""


	gd = GenericDialogPlus("EZ Peeler - Surface segmentation")
	gd.addMessage("""These options define the choices for surface segmentation.
					
					The Minimum threshold slider sets the minimum value at which a voxel is
					recognised as a cell wall.
	
					Interpolation smooths the surface segmentation, then places each vertex an
					equal distance apart. This is useful if segmentation is not smooth. 
					(Recommended off)
	
					Ignoring slices, removes the first N number of XY slices during segmentation,
					removing artefacts introduced by XZ edge detection filters or reflection off
					coverslips. (Recommended off)
				
					**The next step will delete any ROIs in the manager currently**""")
	gd.addCheckbox("""Use default (Otsu method) threshold for each frame? This frame threshold value is """+str(defaultThreshold), True)
	gd.addSlider("Minimum threshold",0, maxPixel, defaultThreshold)
	gd.addSlider("Current visible slice",1, imp.getStackSize(), int(imp.getStackSize()/2))
	gd.addCheckbox("Interpolation", False)
	gd.addNumericField("Interpolation interval ", 8, 0)
	gd.addSlider("Slices to ignore:", 0, int(height),0)

	gd.setModal(False)
	gd.showDialog()
	

	sliders=gd.getSliders()
	
	
	sliderValue=sliders.get(0)
	sliderSlice=sliders.get(1)
	sliceSlid1=sliderSlice.getValue()
	slid1=sliderValue.getValue()

	slid2=slid1
	sliceSlid2=sliceSlid1
	imp.setSlice(sliceSlid1)
	widths, heights= findEpidermis(stack.getPixels(sliceSlid1), width, height, slid1, False, 8, 0, andOp)
	
	#draw a preview segmentation ROI to imp 
	proi = PolygonRoi(widths, heights, len(widths), Roi.POLYLINE)
	proi.setPosition(sliceSlid1)
	proi.setStrokeColor(Color.green)
	overlay=Overlay()
	overlay.add(proi)
	imp.setOverlay(overlay)
	imp.show()
	
	while ((not gd.wasCanceled()) and not (gd.wasOKed())):
		sliders=gd.getSliders()
		sliderValue=sliders.get(0)
		sliderSlice=sliders.get(1)
		sliceSlid1=int(sliderSlice.getValue())
		slid1=sliderValue.getValue()
	
		if (slid2!=slid1 or sliceSlid2!=sliceSlid1):
			slid2=slid1
			sliceSlid2=sliceSlid1
			widths, heights= findEpidermis(stack.getPixels(sliceSlid1), width, height, slid1, False, 8, 0, andOp)
			
			#draw a preview segmentation ROI to imp 
			
			proi = PolygonRoi(widths, heights, len(widths), Roi.POLYLINE)
			proi.setPosition(sliceSlid1)
			proi.setStrokeColor(Color.green)
			overlay=Overlay()
			overlay.add(proi)
			imp.setOverlay(overlay)
			imp.setSlice(sliceSlid1)
			imp.show()
			Thread.sleep(50)
			
			continue
	
	
	
		slid2=slid1
		sliceSlid2=sliceSlid1
		Thread.sleep(50)
	useOtsu = gd.getNextBoolean()
	minThreshold = gd.getNextNumber()
	displaySlice = gd.getNextNumber()
	interpolation = gd.getNextBoolean()
	interpolRes = gd.getNextNumber()
	topSlice = gd.getNextNumber()

	if useOtsu==0:
		minThreshold=defaultThreshold
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 4"
	else:
		canceled = 0
	if gd.wasOKed():
		oked=1
	else:
		oked=0	
		
	return minThreshold, interpolation, canceled, interpolRes, topSlice, oked,useOtsu

def getOptions3(imp, pixels, imp3, epidermisHeightsFull):
	gd = GenericDialogPlus("EZ Peeler -Error correction, epidermis extraction and reslicing")
		
	gd.addMessage(""" You can now check your surface segmentation for errors, and choose any
					autocorrection.
	
					The 'Hole removal' option removes points where the height diverges massively
					from the smoothed surface depth, which can correct segmentation errors

					These options allow the user to specify how deep below the surface to
					peel.
				
					The depth offset allows you start your peel a specified distance after the
					threshold line. This is useful if you don't want noise from the surface
					staining in your stack. The top depth offset allows you to specify a smaller
					offset if your surface goes higher than the top of the stack.
					
					The epidermis thickness allows you to set the thickness of your peel.
					
					Using a Gaussian filter on the render mask, reduces the woodgrain effect of
					images.""")
	gd.addSlider("Current visible slice",1, imp3.getStackSize(), int(imp3.getStackSize()/2))
	gd.addCheckbox("Erode filter peel? (Sometimes buggy with negative values)", False)
	gd.addSlider("y offset (in voxel layers)", -imp3.getHeight(), imp3.getHeight(), 4, 1)
	gd.addNumericField("x and z offset (in voxels for erosion)", 1, 1)
	gd.addSlider("Epidermis y thickness (in voxel layers)",0, imp3.getHeight()*2, 8, 1)
	gd.addNumericField("Epidermis x and z thickness (in voxel layers)", 1, 1)
	gd.addCheckbox("Gaussian filter render mask?", True)
	gd.addCheckbox("Hole removal", True)
	gd.addSlider("Hole divergence threshold", 0, 50, 15, 0.1);
	
	gd.setModal(False)
	gd.showDialog()
	#set up sliders
	
	sliders=gd.getSliders()
	
	sliderSlice=sliders.get(0)
	sliderOffset=sliders.get(1)
	sliderThick=sliders.get(2)
	sliderValue=sliders.get(3)

	sliceSlid1=sliderSlice.getValue()
	offsetSlid1=sliderOffset.getValue()
	thickSlid1=sliderThick.getValue()
	slid1=sliderValue.getValue()
	
	#update slider values
		
	sliceSlid2=sliceSlid1
	offsetSlid2=offsetSlid1
	thickSlid2=thickSlid1
	
	excludedPixels=filter(lambda i: abs(pixels[i]) > slid1/10, xrange(len(pixels)))
	pixelMap=array('i', [0]*len(pixels))
	for i in excludedPixels:
		pixelMap[i]=16674815
	
	#draw a preview of excluded pixels to divergence map 
	excludedProcessor=ColorProcessor(imp.getWidth(), imp.getHeight(), pixelMap)
	proi = ImageRoi(0, 0, excludedProcessor)
	proi.setZeroTransparent(True)
	proi.setOpacity(0.5)
	overlay=Overlay()
	overlay.add(proi)
	imp.setOverlay(overlay)
	imp.show()
	imp.updateAndDraw()

	#update slider values
	slid2=slid1
	sliceSlid2=sliceSlid1
	width=imp.width
	
	for y in xrange(imp.height):
		widths[y]=range(width)
		heights[y]=	epidermisHeightsFull[y*width:(y+1)*(width)]	
	roiYs= [x + offsetSlid1 for x in heights[sliceSlid1-1]] + [x + offsetSlid1 + thickSlid1 for x in reversed(heights[sliceSlid1-1])]
	roiXs= [x for x in widths[sliceSlid1-1]] + [x for x in reversed(widths[sliceSlid1-1])]
	# draw a preview of the segmented slice, of chosen thickness onto imp3 
	proi2 = PolygonRoi(roiXs, roiYs, len(roiXs), Roi.POLYGON) 
	proi2.setPosition(sliceSlid1)
	proi2.setStrokeColor(Color.magenta)

	proi3 = PolygonRoi(xrange(width), heights[sliceSlid1-1], width, Roi.POLYLINE)
	proi3.setPosition(sliceSlid1)
	proi3.setStrokeColor(Color.green)
	
 	overlay2=Overlay()
 	overlay2.add(proi2)
 	overlay2.add(proi3)
 	imp3.setOverlay(overlay2)
 	imp3.show()
	imp3.updateAndDraw()
	imp3.setOverlay(overlay2)
	
	while ((not gd.wasCanceled()) and not (gd.wasOKed())):
		sliders=gd.getSliders()
		sliderSlice=sliders.get(0)
		sliderOffset=sliders.get(1)
		sliderThick=sliders.get(2)
		sliderValue=sliders.get(3)
	
		sliceSlid1=sliderSlice.getValue()
		offsetSlid1=sliderOffset.getValue()
		thickSlid1=sliderThick.getValue()
		slid1=sliderValue.getValue()
		#If slider has changed, update preview	
		if (slid2!=slid1):
			excludedPixels=filter(lambda i: abs(pixels[i]) > slid1/10, xrange(len(pixels)))
			pixelMap=array('i', [0]*len(pixels))
			for i in excludedPixels:
				pixelMap[i]=16674815
			excludedProcessor=ColorProcessor(imp.getWidth(), imp.getHeight(), pixelMap)
			proi = ImageRoi(0, 0, excludedProcessor)
			proi.setZeroTransparent(True)
			proi.setOpacity(0.5)
			overlay=Overlay()
			overlay.add(proi)
			imp.setOverlay(overlay)
			imp.show()
			imp.updateAndDraw()
			slid2=slid1
			
			IJ.selectWindow(imp.getTitle());
			Thread.sleep(50)
			continue
		#If slider has changed, update preview
		if (slid2!=slid1 or sliceSlid2!=sliceSlid1 or offsetSlid2!=offsetSlid1 or thickSlid2!=thickSlid1):
		 	imp3.setSlice(sliceSlid1)
		 	sliceSlid2=sliceSlid1
			offsetSlid2=offsetSlid1
			thickSlid2=thickSlid1
			IJ.selectWindow(imp3.getTitle());
			
			proi2Ys= [x + offsetSlid1 for x in heights[sliceSlid1-1]] + [x + offsetSlid1 + thickSlid1 for x in reversed(heights[sliceSlid1-1])]
			proi2Xs= [x for x in widths[sliceSlid1-1]] + [x for x in reversed(widths[sliceSlid1-1])]
			
			# draw create a new ROI and draw to the mask image image
			proi2 = PolygonRoi(proi2Xs, proi2Ys, len(proi2Xs), Roi.POLYGON) 
			proi2.setPosition(sliceSlid1)
			proi2.setStrokeColor(Color.magenta)

			proi3 = PolygonRoi(xrange(width), heights[sliceSlid1-1], width, Roi.POLYLINE)
			proi3.setPosition(sliceSlid1)
			proi3.setStrokeColor(Color.green)
			
		 	overlay2=Overlay()
		 	overlay2.add(proi2)
		 	overlay2.add(proi3)
		 	imp3.setOverlay(overlay2)
		 	imp3.show()
			imp3.updateAndDraw()
		 	Thread.sleep(50)
		 	
		 	continue
	
	displaySlice = gd.getNextNumber()
	erode = gd.getNextBoolean()
	depthOffset = gd.getNextNumber()
	xzOffset = gd.getNextNumber()
	stackThickness = gd.getNextNumber()
	xzThickness = gd.getNextNumber()
	gaussian = gd.getNextBoolean()
	hdRemoval = gd.getNextBoolean()
	heightDiffMax = gd.getNextNumber()
	if gaussian == 1:
		gaussian = 10.0
	else:
		gaussian = 1
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 4"
	else:
		canceled = 0
	if gd.wasOKed():
		oked=1
	else:
		oked=0	
	return erode, depthOffset, stackThickness, canceled, hdRemoval, heightDiffMax, gaussian, xzOffset, xzThickness, oked

def finalDialog():
	gd = NonBlockingGenericDialog("EZ Peeler - Image check")
		
	gd.addMessage("""Congratulations! your sliced image should now be visible and ready
					for use. If there are errors with your segmentation, hit cancel to
					return to the previous steps and alter your settings accordingly""")
					
	gd.addCheckbox("Keep results image from this time frame?", True)
	gd.addMessage("""If your image was part of a time series, you may now apply the same
	settings to segment the whole time series. Warning: This can take several minutes to process.""")
	gd.addCheckbox("Apply settings to rest of time series?", False)
	gd.addCheckbox("Keep time series 3D stack?", False)
	gd.addCheckbox("Keep time series z projection?", True)
	gd.addCheckbox("Keep time series heightmap?", False)
	gd.addCheckbox("Keep time series areamap?", False)
	gd.showDialog()

	keepPrev = gd.getNextBoolean()
	timeseries = gd.getNextBoolean()
	
	keep3D = gd.getNextBoolean()
	keepZP = gd.getNextBoolean()
	keepHM = gd.getNextBoolean()
	keepAM = gd.getNextBoolean()
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 4"
	else:
		canceled = 0
	if gd.wasOKed():
		oked=1
	else:
		oked=0
	return canceled, timeseries, oked, keep3D, keepZP, keepHM, keepAM, keepPrev	
	
def findEpidermis(pixels, width, height, threshold, interpolation, interpolRes, topSlice,andOp):
	""" iterates through the image, finding the first pixel above the threshold value, then performs user specified error correction"""

	#Create a list of initial values 1 larger than the image height
	
	epidermisHeights=[height+1]*width	  
	
	# Go through pixels, looking for first pixel with a value greater than threshold in the y, placing y coord in a list

	epidermisYs = []
	epidermisXs = []

	
	for x in xrange(width):
		for y in xrange(topSlice, height):
			pixel= width*y+x

			#after finding the first pixel above the threshold value, add the y coord to the list
			if pixels[pixel]& andOp > threshold:
			
				epidermisHeights[x] = y
				#break from looping the y when 1st threshold pixel is found is met -> increases speed drastically! Otherwise need an if statement every loop...
				break


	#Add non-redundant coords to two new lists. Points where no threshold value was met are excluded automatically, to be filled by interpolation later
	#(first and last coords are fixed + therefore added manually, the rest are added non-redundantly in the loop)


	epidermisYs.append(epidermisHeights[0])
	epidermisXs.append(0)

		
	for x in xrange(1, width-1):
		if not epidermisHeights[x] == epidermisHeights[x-1] == epidermisHeights[x+1]:
			if epidermisHeights[x] < height+1:
				epidermisYs.append(epidermisHeights[x])
				epidermisXs.append(x)

		
	epidermisYs.append(epidermisHeights[-1])
	epidermisXs.append(width)

	if epidermisYs[0] > height:
		epidermisYs[0]=epidermisYs[1]
	if epidermisYs[-1] > height:
		epidermisYs[-1]=epidermisYs[-2]		
	
	#return ROI X and Y positions and an array of pure height data for heightmap
	return epidermisXs, epidermisYs

def firstStage(channels, frame, filtertype, smoothsize, canceled1, sobeltype, andOp, GPU):
		#extract the channel you want to base the peeler on
		imp2 = extractChannel(imp1, int(channels)+1, frame, andOp)
		
		#run the chosen filter
		src = smoothFilter(imp2, filtertype, int(smoothsize))
		dst = clij2.create([src.getWidth(), src.getDepth(), src.getHeight()], src.getNativeType())
	
		#reslice the extracted channel and show resultant
		clij2.resliceTop(src, dst)
		src.close()
		imp3 = clij2.pull(dst)
		
		dst.close()
		ip = imp3.getProcessor()
		
		width = ip.getWidth()
		height= ip.getHeight()
	
		# Redraw based on max pixel value
		findEdge(imp3, sobeltype)
		stats =StackStatistics(imp3) 
		frame1Mean=stats.mean
		
		IJ.setMinAndMax(imp3, 0, stats.max)
		IJ.run(imp3, "Grays", "stack")
	
		imp3.setTitle("Resliced channnel "+ str(channels)+ ", frame "+ str(frame))
		imp3.show()
		
		IJ.setAutoThreshold(imp3, "Otsu dark stack")
		defaultThreshold=imp3.getProcessor().getMinThreshold()
		IJ.resetThreshold(imp3)
		stack = imp3.getStack()
		return imp2, imp3, width, height, stats, stack, defaultThreshold
def secondStage(minThreshold, interpolation, canceled2,  interpolRes, topSlice, useOtsu, defaultThreshold):
	if useOtsu==1:
			minThreshold=defaultThreshold
	
	heights = {}
	widths = {}
	


	epidermisHeightsFull = []
	xvertices= [range(width)]
	
	for sliceN in xrange(imp3.getNSlices()):


		#get x + y positions of the top of the epidermis for this slice
		widths[sliceN+1], heights[sliceN+1] = findEpidermis(stack.getPixels(sliceN+1), width, height, minThreshold, interpolation, interpolRes, topSlice, andOp)

		
		#add to ROI and editor
		proi = PolygonRoi(widths[sliceN+1], heights[sliceN+1], len(widths[sliceN+1]), Roi.POLYLINE) 
		
		if interpolation:
			proi2 = proi.getInterpolatedPolygon(interpolRes, 1)
		
			n_points = proi2.npoints
			x = proi2.xpoints
			y = proi2.ypoints
			x[0]=0
			x[-1] = width
			proi = PolygonRoi(x,y,n_points, Roi.POLYLINE)
			
		proi.setPosition(sliceN+1)
		#By extracting the contained pixels in the multipoint line roi, we interpolates all the missing data in the heightmap.
		points = proi.getContainedPoints()
		points= sorted(points, key=lambda Point:Point.x)
		yvertices=[points[0].y]
		
		#create new, non redundant list of heights at 1 space interval
		for i in xrange(1, len(points)):
			if points[i].x != points[i-1].x:
				yvertices.append(points[i].y)
		while len(yvertices) > width:
			yvertices.pop(-1)
		
		epidermisHeightsFull+=yvertices
	print epidermisHeightsFull
	
	print "points"
	print points
	
	
	#make heightmap
	heightMapArray = array( "f", epidermisHeightsFull)
	fp= FloatProcessor(width, len(heightMapArray)/width, epidermisHeightsFull, None)
	heightsImp= ImagePlus("Uncalibrated Heightmap", fp)
	heightsImp.show()
	
	fp2= FloatProcessor(width, len(heightMapArray)/width, epidermisHeightsFull, None)
	
	#Create a height 'Divergence map' and use this to remove peaks/holes in the segmentation

	GaussianBlur().blurGaussian(fp2, 30)
	blurredHeights= ImagePlus("Blurred heightmap " + str(frame), fp2)
	imgBlur=ImgLib.wrap(blurredHeights)

	imgHeights=ImgLib.wrap(heightsImp)
	sub = Compute.inFloats(Subtract(imgBlur, imgHeights))
	subIP= ImgLib.wrap(sub).getProcessor()
	subPixels=subIP.getPixelsCopy()
	#blurredHeights.show()
	impSub= ImgLib.wrap(sub)
	impSub.setTitle("Divergence map")
	IJ.setMinAndMax(impSub, min(subPixels), max(subPixels))	
	return heights, widths, epidermisHeightsFull, xvertices, heightMapArray, fp, heightsImp, blurredHeights, imgBlur, imgHeights, sub, subIP, subPixels, impSub


def thirdStage(epidermisHeightsFull, erode, depthOffset, stackThickness, canceled3, hdRemoval, heightDiffMax, gaussian,xzOffset, xzThickness):

	widths={}
	heights={}
	if hdRemoval:
		
	#Create a dict of lists, for each list, copy heights from heightimage, leaving out choords that exceed the threshold in the divergence map.
		excludedPixels=filter(lambda i: abs(subPixels[i]) > heightDiffMax, xrange(len(subPixels)))
		excludeMapArray=array('b', [0]*len(subPixels))
		for i in excludedPixels:
			excludeMapArray[i]=100
		excludeP= ByteProcessor(width, len(subPixels)/width, excludeMapArray)
		excludeP.erode()
		excludeP.erode()
		excludeP.setThreshold(90, 110, ImageProcessor.NO_LUT_UPDATE)
		excludedImp=ImagePlus("exclusion map", excludeP)
		excludeROI = ThresholdToSelection.run(excludedImp)
		excludedImp.close()
		correctedHeightsImp = ImagePlus("corrected Heightmap", fp).duplicate()
		correctedHeightsImp.setRoi(excludeROI)
		
		#blur in missing values with gaussian
		IJ.run(correctedHeightsImp, "Gaussian Blur...", "sigma=50")
		for x in range(10):
			IJ.run(correctedHeightsImp, "Gaussian Blur...", "sigma=2")
		for x in range(10):
			IJ.run(correctedHeightsImp, "Gaussian Blur...", "sigma=1.5")
		correctedHeightsP=correctedHeightsImp.getProcessor()
		epidermisHeightsFull=correctedHeightsP.getPixels()
		correctedHeightsImp.close()
	for y in xrange(impSub.height):
		widths[y]=range(width)
		heights[y]=	epidermisHeightsFull[y*width:(y+1)*(width)]

	


	#create an image to draw the ROIs to
	binaryMask = IJ.createImage("Binary Mask", "8-bit black", width, height,(len(heights)))
	
	stack2=binaryMask.getImageStack()
	
	#iterate through the ROI chords, and draw to the new image
	if erode:
		for sliceN in range(len(heights)):
			# make new ROIs, the size and thickness requested by the user	
			roiYs= [x+1 for x in heights[sliceN]] +[height, height]
			roiXs= [x for x in widths[sliceN]] + [width+1, 0]

			# create a new ROI and draw to the mask image image
			proi = PolygonRoi(roiXs, roiYs, len(roiXs), Roi.POLYGON) 
			proi.setPosition(sliceN+1)

			bsip = stack2.getProcessor(sliceN+1)
			bsip.setColor(gaussian)
			bsip.fill(proi)
			

		BS2 = ImagePlus("Binary Stack", stack2)
		IJ.setMinAndMax(BS2, 0, 1)
		
		src2 = clij2.push(BS2)
		BS2.close()	
		dst2 = clij2.create(src2)

		clij2.minimum3DSphere(src2, dst2, int(xzOffset),int(depthOffset),int(xzOffset))
		clij2.minimum3DSphere(dst2, src2, int(xzThickness), int(stackThickness), int(xzThickness))

		dst3=clij2.create(src2)

		clij2.subtract(dst2, src2, dst3)

		BS = clij2.pull(dst3)
		#BS.show()
		src2.close()
		dst2.close()
		dst3.close()
	else:
		for sliceN in range(len(heights)):
			# make new ROIs, the size and thickness requested by the user	
			roiYs= [x + depthOffset for x in heights[sliceN]] + [x + depthOffset + stackThickness for x in reversed(heights[sliceN])]
			roiXs= [x for x in widths[sliceN]] + [x for x in reversed(widths[sliceN])]

			# draw create a new ROI and draw to the mask image image
			proi = PolygonRoi(roiXs, roiYs, len(roiXs), Roi.POLYGON) 
		
			proi.setPosition(sliceN+1)
		
			bsip = stack2.getProcessor(sliceN+1)
		
			bsip.setColor(gaussian)
		
			bsip.fill(proi)

		BS = ImagePlus("Binary Stack", stack2)




	#reslice the mask image, so that it can be applied to the original image 
	
	binaryMask=Slicer().reslice(BS)
	stack2 = binaryMask.getImageStack()
	stack3= ImageStack(imp1.width, imp1.height)
	#if wanted, apply a blur to the mask, to prevent woodgrain aliasing errors

	for i in xrange(imp1.getNSlices()):
		channelI = stack2.getProcessor(i+1)
		for j in xrange((imp1.getNChannels())):
			stack3.addSlice(None, channelI)	
			
	heightmapFinal=[]
	

	#make heightmap
	cal = imp1.getCalibration()  
	depth= cal.pixelDepth
	heightmapFinal = map(lambda x:x*depth, epidermisHeightsFull)
	heightMapArray = array( "f", heightmapFinal)
	fp3= FloatProcessor(width, len(heightMapArray)/width, heightMapArray, None)
	heightsImp2= ImagePlus("Calibrated heightmap", fp3)
	heightsImp2.setCalibration(cal)
	heightsImp2.show()
	areaImp = areamap( heightsImp2, cal.pixelWidth, cal.pixelHeight)
	areaImp.setTitle("area map")
	areaImp.setCalibration(cal)
	areaStat =areaImp.getStatistics()
	#IJ.setMinAndMax(areaImp, areaStat.min, areaStat.max)
	
	areaImp.show()
	try:
		IJ.run(areaImp, "16_colors", "")
	except: print "bugger"
	imp6 = ImagePlus("Binary mask to channels", stack3)
	imp5 = extractFrame(imp1, frame)
	imp5.setCalibration(cal)
	imp6.setCalibration(cal)
	if gaussian == 10.0:
		IJ.run(imp6, "Gaussian Blur...", "sigma=1 stack")
		imp4 = ImageCalculator().run("Multiply create 32-bit stack" , imp5, imp6)
		IJ.run(imp4, "Divide...", "value=10.0 stack");
	else:
		imp4 = ImageCalculator().run("Multiply create stack" , imp5, imp6)
	
	
	stats =StackStatistics(imp4) 
	IJ.setMinAndMax(imp4, 0, stats.max)
	IJ.run(BS, "Grays", "stack")
	IJ.setMinAndMax(BS, 0, 1)
	IJ.setMinAndMax(binaryMask, 0, 1)
	imp4.setTitle("Segmented surface")
	imp4.show()
	imp5.close()
	
	sumProjImp = ZProjector.run(imp4, "sum")
	sumProjImp.show()
	return widths, heights, epidermisHeightsFull, stack3, heightsImp2, imp4, imp6, areaImp, sumProjImp

"""********************************actual script*****************************************"""



#get the current image
imp1= IJ.getImage()
#ImageConverter(imp1).convertToGray16() 
finish =0

canceled1 =0
canceled2 =0
canceled3 =0
canceled4 =0

oked1 = 0
oked2 = 0
oked3 = 0
oked4 = 0
stage=1

headless = 0
keep3D=1
keepZP=1
keepHM=1
keepAM=1
while (canceled1 == 0 and oked4==0):
	oked2 = 0
	oked3 = 0
	oked4 = 0
	canceled2 =0
	canceled3 =0
	canceled4 =0
	while stage==1:
		#cleanup open images
		try:
			imp2.close()
		except:
			print "imp2 already closed"
		try:
			imp3.close()
		except:
			print "imp3 already closed"
		
		#ask the user the settings they want
		userOptions1 = getOptions1(imp1)
	
		channels, frame, filtertype, smoothsize, canceled1, sobeltype, andOp, GPU, oked1 = userOptions1
		if canceled1==1:
			stage=0
			break
		if oked1==1:
			stage=2
		
		frame=int(frame)
		clij2.getInstance(GPU)
		

		imp2, imp3, width, height, stats, stack, defaultThreshold = firstStage(channels, frame, filtertype, smoothsize, canceled1, sobeltype, andOp, GPU)


		
	#run the second dialog
	while stage==2:
		oked3 = 0
		oked4 = 0
		canceled3 =0
		canceled4 =0
		#cleanup open images
		try:
			heightsImp.close()
		except:
			print "heightsImp already closed"
		try:
			impSub.close()
		except:
			print "impSub already closed"			
		
		userOptions2 = getOptions2(imp3, stats.max, height-1, stack, andOp, width, defaultThreshold)	
		
			
		#minThreshold, interpolation, hdRemoval, heightDiffMax, canceled2,  interpolRes, topSlice = userOptions2
		minThreshold, interpolation, canceled2,  interpolRes, topSlice, oked2 , useOtsu= userOptions2

		
		if canceled2==1:
			stage=1
			break
		if oked2==1:
			stage=3
		heights, widths, epidermisHeightsFull, xvertices, heightMapArray, fp, heightsImp, blurredHeights, imgBlur, imgHeights, sub, subIP, subPixels, impSub=secondStage(minThreshold, interpolation, canceled2,  interpolRes, topSlice, useOtsu, defaultThreshold)

	while stage==3:

		try:
			BS.close()
		except:
			print "BS already closed"
		try:
			BS2.close()
		except:
			print "BS2 already closed"
		try:
			imp4.close()
		except:
			print "imp4 already closed"
		try:
			imp5.close()
		except:
			print "imp5 already closed"
		try:
			imp6.close()
		except:
			print "imp6 already closed"
		try:
			areaImp.close()				
		except:
			print "areaImp already closed"
			
		try:
			impSub.close()				
		except:
			print "impSub already closed"
		try:
			heightsImp2.close()
		except:
			print "heightsImp2 already closed"
		oked4 = 0
		canceled4 =0
		

		
		userOptions3 = getOptions3(impSub, subPixels, imp3, epidermisHeightsFull)

		erode, depthOffset, stackThickness, canceled3, hdRemoval, heightDiffMax, gaussian,xzOffset, xzThickness, oked3= userOptions3
		
		
		if canceled3==1:
			stage=2
			break
		if oked3==1:
			stage=4

		widths, heights, epidermisHeightsFull, stack3, heightsImp2, imp4, imp6, areaImp, sumProjImp = thirdStage(epidermisHeightsFull, erode, depthOffset, stackThickness, canceled3, hdRemoval, heightDiffMax, gaussian,xzOffset, xzThickness)	

		

		
	#Final dialog loop
		
	while stage==4:
		finalOptions=finalDialog()
		canceled4, timeseries, oked4, keep3D, keepZP, keepHM, keepAM, keepPrev = finalOptions
		if canceled4==1:
			stage=3
			break
		if oked4==1:
			stage=5

imp2.close()
imp3.close()
impSub.close()
heightsImp.close()

if keepPrev == 0:
	sumProjImp.close()
	heightsImp2.close()
	imp4.close()
	areaImp.close()
	imp6.close()
if timeseries == 1:
	headless = 0
	stack4= ImageStack(imp1.width, imp1.height)
	stack5 = ImageStack(imp1.width, imp1.height)
	stack6 = ImageStack(imp1.width, imp1.height)
	stack7 = ImageStack(imp1.width, imp1.height)

	for frame in xrange(1, imp1.getNFrames()+1):
		print "Frame:  " + str(frame)
		imp2, imp3, width, height, stats, stack, defaultThreshold = firstStage(channels, frame, filtertype, smoothsize, canceled1, sobeltype, andOp, GPU)
		heights, widths,  epidermisHeightsFull, xvertices, heightMapArray, fp, heightsImp, blurredHeights, imgBlur, imgHeights, sub, subIP, subPixels, impSub=secondStage(minThreshold, interpolation, canceled2,  interpolRes, topSlice, useOtsu, defaultThreshold)
		widths, heights, epidermisHeightsFull, stack3, heightsImp2, imp4, imp6, areaImp, sumProjImp = thirdStage(epidermisHeightsFull, erode, depthOffset, stackThickness, canceled3, hdRemoval, heightDiffMax, gaussian,xzOffset, xzThickness)	
		if keepAM:
			stack7.addSlice(areaImp.getProcessor())
		if keepHM:
			stack6.addSlice(heightsImp2.getProcessor())
		if keep3D:
			imp4stack=imp4.getImageStack()
			for i in xrange(1, imp4stack.getSize()+1):	
				try:	
					stack4.addSlice(imp4stack.getProcessor(i))	
				except: print "FAIL"
		if keepZP:
			imp5stack = sumProjImp.getImageStack()
			for i in xrange(1, imp5stack.getSize()+1):	
				try:	
					stack5.addSlice(imp5stack.getProcessor(i))	
				except: print "FAIL"
		
		sumProjImp.close()		
		imp6.close()
		areaImp.close()
		heightsImp2.close()
		imp4.close()
		heightsImp.close()	
		impSub.close()
		imp2.close()
		imp3.close()
		
	if keep3D:
		imp7 = ImagePlus("Segmented timeseries", stack4)
		imp7.setDimensions(imp1.getNChannels(), imp1.getNSlices(), imp1.getNFrames())
		imp7 = CompositeImage(imp7, CompositeImage.COMPOSITE)  
		imp7.show()
	if keepAM:
		imp10 = ImagePlus("Areamap timeseries", stack7)
		imp10.setDimensions(1, 1, imp1.getNFrames())
		imp10.show()
	
	if keepHM:
		imp9 = ImagePlus("Segmented timeseries heightmaps", stack6)
		imp9.setDimensions(1, 1, imp1.getNFrames())
		imp9.show()
	if keepZP:
		imp8 = ImagePlus("Segmented timeseries sum projection", stack5)
		imp8.setDimensions(imp1.getNChannels(), 1, imp1.getNFrames())
		imp8 = CompositeImage(imp8, CompositeImage.COMPOSITE)  
		imp8.show()

