"""
*************************************************************************************************
							Epidermal Z stack Peeler (EZ Peeler)
											v 0.16

						Written by Jim Rowe (Stuart Casson's lab, Sheffield, UK)
		In case it's not obvious, I'm a biologist, not a coder, so apologies for sloppy spaghetti code
						
										@BotanicalJim
									james.rowe@slcu.cam.ac.uk
								orcid.org/0000-0002-3523-4347

	Development started 15/05/2017 when I should have been doing physiology work.... 

**************************************************************************************************
"""


#import libraries required for plugin

from ij import IJ, ImagePlus, ImageStack
from ij.gui import DialogListener, Roi, PolygonRoi, NonBlockingGenericDialog
from ij.process import ImageProcessor, StackStatistics, ImageConverter
from ij.plugin import Slicer, ChannelSplitter, ImageCalculator 
from ij.plugin.filter import GaussianBlur, RankFilters
from ij.plugin.frame import RoiManager

"""**********************************define functions*******************************************"""

def extractChannel(imp, nChannel, nFrame):
	"""extract a channel from the image, returning a new 16 bit imagePlus labelled with the channel name"""
	
	imp2=ChannelSplitter.getChannel(imp, nChannel)

	imp3 = ImagePlus("Channel " + str(nChannel), imp2).duplicate()
	
	ImageConverter(imp3).convertToGray16() 
	return imp3




def smoothFilter(imp, filtertype, radius):
	"""iterate through slices to perform either gaussian or median filter"""
	if filtertype == "Gaussian":
		for slice in xrange(imp.getStackSize()):
			GaussianBlur().blurGaussian(imp.getStack().getProcessor(slice+1), radius, radius, 0.01)
	elif filtertype == "Median":
		for slice in xrange(imp.getStackSize()):
			RankFilters().rank(imp.getStack().getProcessor(slice+1), radius, RankFilters.MEDIAN)
	#return as imagePlus
	return imp



	
def getOptions1(imp):
	"""Get user defined options for image preprocessing"""

	gd = NonBlockingGenericDialog("EZ Peeler - Preprocessing")

	gd.addMessage("""                     Written by Jim Rowe
	
					These options define what preprocessing you want to apply 
					to the image for the purpose of segmenting the epidermis.
		
					These filters will not be applied to the final output and
					are only for the purpose of segmentation. Smoothing values
					work well at removing granular noise from confocal stacks.
					Edge detection uses a Sobel filter to improve surface 
					contrast between the surface and background""")

	
	#create a list of the channels in the provided imagePlus
	types = []
	for i in xrange(imp.getNChannels()):
		types.append(str(i))
	#user can pick which channel to base the peeling on
	
	gd.addChoice("Channel to use for segmentation", types, types[-1])


	#user specifies the smoothing process and the size of the smoothing
	smoothtype = ["Gaussian", "Median", "none"]
	gd.addChoice("Image smoothing type", smoothtype, smoothtype[0])
	
	smoothing=["1", "2", "3", "4", "5", "6"]
	gd.addChoice("Image smoothing sigma", smoothing, smoothing[3])
	#Use specifies whether they want to use edge detection
	sobel=["none", "Sobel"]
	gd.addChoice("Edge detection filter", sobel, sobel[1])
	 
	gd.showDialog()
	

	channel = gd.getNextChoice()
	filtertype = gd.getNextChoice()
	smoothsize = gd.getNextChoice()
	sobeltype = gd.getNextChoice()
	
	# if dialog canceled then return nothing
	if gd.wasCanceled():
		finish = 1
		print "canceled dialog 1"
	else:
		finish = 0	
	return channel, filtertype, smoothsize, finish, sobeltype


def getOptions2(imp, maxPixel):

	"""Get user defined options for segmentation"""


	gd = NonBlockingGenericDialog("EZ Peeler - Surface segmentation")
	gd.addMessage("""These options define the choices for surface segmentation.
					
					The Minimum threshold slider sets the minimum value at which a voxel is
					recognised as a cell wall.
	
					Interpolation smooths the surface segmentation, then places each vertex an
					equal distance apart. This is useful of segmentation requires a lot of manual
					error correction, otherwise (Recommended off)
	
					The 'Hole removal' option removes points where the height with a sudden
					drop in surface depth, which can correct errors where the epidermal surface
					is above the top of the stack, or the surface at that coordinate hasn't hit
					the threshold set. (Recommended on)
				
					**The next step will delete any ROIs in the manager currently**""")

	gd.addSlider("Minimum threshold",0, maxPixel, int(maxPixel*0.1))
	gd.addCheckbox("Interpolation", False)
	gd.addNumericField("Interpolation interval ", 8, 0)
	gd.addCheckbox("Hole removal", True)
	gd.addNumericField("Maximum allowable hole depth (in voxels) ", 4, 0)


	
	gd.showDialog()
	
	minThreshold = gd.getNextNumber()
	interpolation = gd.getNextBoolean()
	interpolRes = gd.getNextNumber()
	hdRemoval = gd.getNextBoolean()
	heightDiffMax = gd.getNextNumber()

	
	if gd.wasCanceled():
		canceled1 = 1
		print "canceled dialog 2"
	else:
		canceled1 = 0	
		
	return minThreshold, interpolation, hdRemoval, heightDiffMax, canceled1, interpolRes

def getOptions3():
	gd = NonBlockingGenericDialog("EZ Peeler - Epidermis extraction and reslicing")
		
	gd.addMessage(""" You can now check your surface segmentation for errors, and correct them
					by using the ROI manager.
	
					These options allow the user to specify how deep below the surface to
					peel.
				
					The depth offset allows you start your peel a specified distance after the
					threshold line. This is useful if you don't want noise from the surface
					staining in your stack. The top depth offset allows you to specify a smaller
					offset if your surface goes higher than the top of the stack.
					
					The epidermis thickness allows you to set the thickness of your peel.
					
					Using a Gaussian filter on the render mask, reduces the woodgrain effect of
					images.""")
					
	gd.addNumericField("Depth offset (in voxel layers)", 5, 0)
	gd.addNumericField("Top depth offset (in voxel layers)", 3, 0)
	gd.addNumericField("Epidermis thickness (in voxel layers)", 12, 0)
	gd.addCheckbox("Gaussian filter render mask?", True)

	
	gd.showDialog()
	
	depthOffset = gd.getNextNumber()
	topOffset = gd.getNextNumber()
	stackThickness = gd.getNextNumber()
	gaussian = gd.getNextBoolean()
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 3"
	else:
		canceled = 0	
	return depthOffset, stackThickness, canceled, gaussian, topOffset

def finalDialog():
	gd = NonBlockingGenericDialog("EZ Peeler - Image check")
		
	gd.addMessage("""Congratulations! your sliced image should now be visible and ready
					for use. If there are errors with your segmentation, hit cancel to
					return to the previous steps and alter your settings accordingly""")
	
	gd.showDialog()
	
	if gd.wasCanceled():
		canceled = 1
		print "canceled dialog 4"
	else:
		canceled = 0
	return canceled
	
def findEpidermis(pixels, width, height, threshold, interpolation, interpolRes, hdRemoval, heightDiffMax):
	""" iterates through the image, finding the first pixel above the threshold value, then performs user specified error correction"""

	#Create a list of initial values 1 larger than the image height
	
	epidermisHeights=[height+1]*width	  
	
	# Go through pixels, looking for first pixel with a value greater than threshold in the y, placing y coord in a list

	epidermisYs = []
	epidermisXs = []

	for x in xrange(width):
		for y in xrange(height):
			pixel= width*y+x

			#after finding the first pixel above the threshold value, add the y coord to the list
			if pixels[pixel]>threshold:
			
				epidermisHeights[x] = y
					
				#break from looping the y when 1st threshold pixel is found is met -> increases speed drastically! Otherwise need an if statement every loop...
				break

	#set any points where no threshold was reached to be the same as their neighbour (auto-removes errors)
	for x in range(2,len(epidermisHeights)-1):
		if epidermisHeights[x] >= height+1:
			epidermisHeights[x]=epidermisHeights[x-1]
	
	if  hdRemoval== True:
		for x in xrange(2, len(epidermisHeights)):
			if epidermisHeights[x]-epidermisHeights[x-1] > heightDiffMax:
					if epidermisHeights[x-1] < height+1:
						epidermisHeights[x]=epidermisHeights[x-1]
						
						
		for x in xrange((len(epidermisHeights)-2), 0, -1):
			if epidermisHeights[x]-epidermisHeights[x+1] > heightDiffMax:
					if epidermisHeights[x+1] < height+1:
						epidermisHeights[x]=epidermisHeights[x+1]
						
	
	

	#Add non-redundant coords to two new lists
	#(first and last coords are fixed + therefore added manually, the rest are added non-redundantly in the loop)


	epidermisYs.append(epidermisHeights[0])
	epidermisXs.append(0)

		
	for x in xrange(1, width-1):
		if not epidermisHeights[x] == epidermisHeights[x-1] == epidermisHeights[x+1]:
			epidermisYs.append(epidermisHeights[x])
			epidermisXs.append(x)
		
	epidermisYs.append(epidermisHeights[-1])
	epidermisXs.append(width)

	
	return epidermisXs, epidermisYs



	
"""********************************actual script*****************************************"""


#get the current image
imp1= IJ.getImage()
finish =0

while finish == 0:
	canceled2=0
	canceled3=0
	#ask the user the settings they want
	userOptions1 = getOptions1(imp1)

	channels, filtertype, smoothsize, finish, sobeltype = userOptions1
	#extract the channel you want to base the peeler on
	#if the user doesn't cancel run the script

	if finish ==1:
		break
	
	imp2 = extractChannel(imp1, int(channels)+1, 1)
	#run the chosen filter
	smoothFilter(imp2, filtertype, int(smoothsize))

	imp2.show()

	#reslice the extracted channel and show resultant
	imp3=Slicer().reslice(imp2)


	ip = imp3.getProcessor()
	
	width = ip.getWidth()
	height= ip.getHeight()	

	# Redraw based on max pixel value
	stats =StackStatistics(imp2) 
	IJ.setMinAndMax(imp3, 0, stats.max)
	IJ.run(imp3, "Grays", "stack")
	if sobeltype =="Sobel":
		IJ.run(imp3, "Find Edges", "")
	imp3.show()
		
	stack = imp3.getStack()	
	#run the second dialog
	while canceled2 ==0:
			canceled3 = 0
			userOptions2 = getOptions2(imp3, stats.max)	
			minThreshold, interpolation, hdRemoval, heightDiffMax, canceled2,  interpolRes = userOptions2
			
			if canceled2 ==1:
				break
			#open the ROI mangager and delete any existing ROIs
			
			heights = {}
			widths = {}
			rm = RoiManager.getInstance()
		
			if not rm:
				rm = RoiManager()	
	
			rm.runCommand("Associate", "true")
			rm.runCommand("Deselect")
			rm.runCommand("Delete") 
		
			for sliceN in xrange(imp3.getNSlices()):

				#get x + y positions of the top of the epidermis for this slice
				widths[sliceN+1], heights[sliceN+1] = findEpidermis(stack.getPixels(sliceN+1), width, height, minThreshold, interpolation, interpolRes, hdRemoval, heightDiffMax)
				#add to ROI and editor
				proi = PolygonRoi(widths[sliceN+1], heights[sliceN+1], len(widths[sliceN+1]), Roi.POLYLINE) 
				if interpolation:
					proi2=proi.getInterpolatedPolygon(interpolRes, 1)
					proi = PolygonRoi(proi2, Roi.POLYLINE)
				proi.setPosition(sliceN+1)
				rm.addRoi(proi)
			
			#show all active ROIs in the manager
			rm.runCommand(imp3, "Show All")	

			#run dialog 3
			
			while canceled3 is not 1:
			
				userOptions3 = getOptions3()
		
				depthOffset, stackThickness, canceled3, gaussian, topOffset = userOptions3
				offsetDiff = topOffset - depthOffset
				if canceled3 ==1:
					break
				
				#copy ROI coords from manager to dicts of lists
				rois = rm.getInstance().getRoisAsArray()
				widths={}
				heights ={}
				j=0
				for i in rois:
					#get the X and Y coords
					widths[j] = i.getXCoordinates()
					heights[j] = i.getYCoordinates()
					#add on the offset to the X and Y coords, to give the true vertex position
					widths[j]= map(lambda x:x+i.getXBase(), widths[j])
					heights[j]= map(lambda x:x+i.getYBase(), heights[j])
					j = j + 1

				#create an image to draw the ROIs to
				binaryMask = IJ.createImage("Binary Mask", "32-bit black", width, height,(len(heights)))
				stack2=binaryMask.getImageStack()

				#iterate through the ROI chords, and draw to the new image
				for sliceN in range(len(heights)):
					# make new ROIs, the size and thickness requested by the user	
					roiYs= [x + depthOffset for x in heights[sliceN]] + [x + depthOffset + stackThickness for x in reversed(heights[sliceN])]
					"""
					for y in range(len(roiYs)/2):
						if roiYs[y] == depthOffset:
							roiYs[y]= roiYs[y] + offsetDiff
					
					for y in range(len(roiYs)/2, len(roiYs)):
						if roiYs[y] == depthOffset + stackThickness:
							roiYs[y] = roiYs[y]  + offsetDiff
					"""
					roiXs= [x for x in widths[sliceN]] + [x for x in reversed(widths[sliceN])]

					# draw create a new ROI and draw to the mask image image
					proi = PolygonRoi(roiXs, roiYs, len(roiXs), Roi.POLYGON) 
				
					proi.setPosition(sliceN+1)
				
					bsip = stack2.getProcessor(sliceN+1)
				
					bsip.setColor(10.0)
				
					bsip.fill(proi)

				BS = ImagePlus("Binary Stack", stack2)

				#reslice the mask image, so that it can be applied to the original image 
				
				binaryMask=Slicer().reslice(BS)
				stack2 = binaryMask.getImageStack()
				stack3= ImageStack(imp1.width, imp1.height)
				#if wanted, apply a blur to the mask, to prevent woodgrain aliasing errors
				
				

				

				for i in xrange(imp1.getNSlices()):
					if gaussian:
						GaussianBlur().blurGaussian(stack2.getProcessor(i+1), 2,2, 0.01)
					channelI = stack2.getProcessor(i+1)
					for j in xrange((imp1.getNChannels())):
						stack3.addSlice(None, channelI)	
	
				imp3 = ImagePlus("Binary mask to channels", stack3)
				imp4 = ImageCalculator().run("Multiply create stack" , imp1, imp3)
				stats =StackStatistics(imp4) 
				IJ.setMinAndMax(imp4, 0, stats.max)
				IJ.run(BS, "Grays", "stack")
				BS.show()
				binaryMask.show()
				imp4.show()


				#dialog decisions
				
				canceled4 = 0

				while canceled4 == 0 :
					canceled4 = finalDialog()
					
					if canceled4 == 1:
						break
					canceled3 = 1
					finish = 1
					canceled2 = 1
					break

print "Plugin finished"



