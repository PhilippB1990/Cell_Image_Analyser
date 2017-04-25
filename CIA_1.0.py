########################################################
#### Cell-Image-Analyser (CIA)		 				####
####												####
#### created by: 	Philipp Borchert				####
#### contact:		pborchert.dev@gmail.com			####
#### Used software: Fiji / Imagej2 and its API		####
####												####
########################################################


##########################
##### import classes #####
##########################

# ImageStack, CompositeImage for RGB-Stack
from ij import IJ, ImagePlus,ImageStack, CompositeImage
from ij.process import ImageStatistics as IS
import os
from ij.io import DirectoryChooser

#necessary for inverting (and many more things :p)
from ij.process import ImageProcessor,ImageConverter
#to enhance contrast
from ij.plugin import ContrastEnhancer as CE
#... well, ImageCalculation :p
from ij.plugin import ImageCalculator as IC
#
from ij.process import FloatProcessor
# for watershed + remove outliers (RankFilters)
from ij.plugin.filter import EDM, ParticleAnalyzer as PA,RankFilters as RF
#resultsTable for output
from ij.measure import ResultsTable
#roiManager
from ij.plugin.frame import RoiManager
#to check calculation time
from datetime import datetime
#for math.fabs
import math
#to get an dialog menu (Interface)
from ij.gui import GenericDialog  
#import to get access to "infinity" with Double.POSITIVE_INFINITY
from java.lang import Double
#import math functions for every kind of calculation
from math import sqrt	# square root
from math import pow	# power
from math import fmod	# modulo
from math import floor	# floor (to lower value)
#import sorting option with itemgetter
from operator import itemgetter


options = IS.MEAN | IS.MEDIAN | IS.MIN_MAX | IS.AREA

#global variables // use "global var" in function to link to global variables
##################################
#_______________________________________________________ global variables _______________________________________________________

# save "mean"-values in dictionary
meanValues 			= {}
# save count-values with CountDic
CountDic 			= {}
# needed in partcile analyzer (global)
StatTable 			= ResultsTable()

#Training variables
xml_list_X 			= []
xml_list_Y 			= []
PA_list_X 			= []
PA_list_Y			= []
PA_list_Distance	= []
#store start-values for output
start_value_array 	= []
#manual_auto_diff	= 0
manual_hitlist		= []
#exception variables
broken_dict 		= {}
#error string (which errors occured during processing)
error_dict 			= {}

#initiate error_dict
#file input errors (FIE)
error_dict['FIE'] 	= []
#parameter input error (PIE)
error_dict['PIE'] 	= []
#output error (OE)
error_dict['OE'] 	= []
#Iteration_ThresHold values
orgTH 				= 0
newTH 				= 0
storeTH				= 0
#saves distance to 255 (for adjusted th)
distance_TH			= 0
#iteration step size
iteration_step_size	= 0
#adjust value for TH
adjust 				= 0
#
oldAvPartSize   	= 0


#script versions
pubVersion 			= 1.0
Feature_Version 	= 1.0
CellCounterVersion 	= 2.2

#_______________________________________________________ Method- / Function- area ___________________________________________________________

# Website fore Imagej API: http://rsb.info.nih.gov/ij/developer/api/overview-summary.html

###############################
##### statistics function #####
###############################
# statistics : http://rsb.info.nih.gov/ij/developer/api/ij/process/ImageStatistics.html

def getStatistics(imp):
	#return statistics for the given ImagePlus object
	global options
	ip = imp.getProcessor()
	stats = IS.getStatistics(ip, options, imp.getCalibration())
	return stats.mean
	

###########################
##### invert function #####
###########################
# imageprocessor :http://rsb.info.nih.gov/ij/developer/api/ij/process/ImageProcessor.html

def invertImage(image):

	ip = image.getProcessor()
	ip.invert()
	image = ImagePlus("inv", ip)


	return image


#################################################
##### image into stacks into image function #####
#################################################			

#gets image and int value to return the right value
def getStacks(image,choice):
	stack = image.getImageStack()
	for i in xrange(1,image.getNSlices()+1):
		#get slice at index i
		colourProcessor		= stack.getProcessor(i)
		#colours
		red 	= colourProcessor.toFloat(0,None)
		green 	= colourProcessor.toFloat(1,None)
		blue 	= colourProcessor.toFloat(2,None)
	if 		choice == 1:
		return red, green, blue
	elif 	choice == 2:
		return green, blue
	elif	choice == 3:
		return red, green
	elif 	choice == 4:
		return red
	else:
		return green
		

#################################
##### Image Filter function #####
#################################
# use this function to alternate your image if you got too much noise (Might be worth adding gausian blure or similar)

def filterImage(image,x):
	ip = image.getProcessor()
	for i in xrange(0,x+1):
		ip.smooth()
	return ip


	
#######################################
##### Watershed on Image function #####
#######################################
	
def watershed(image):
	ip = image.getProcessor()
	EDM().toWatershed(ip)

#####################################
##### ParticleAnalyzer functions #####
#####################################

# source: http://rsb.info.nih.gov/ij/developer/api/ij/plugin/filter/ParticleAnalyzer.html

# use particle analyser to count particle in image

def ParAna(image,table,CountDictionary,filename,imp):


	options = PA.SHOW_ROI_MASKS \
		+ PA.INCLUDE_HOLES  \
		+ PA.SHOW_OUTLINES \

	roim = RoiManager(True)
	
	#added PA.CENTROID to get centroids of all particles (X,Y Values)
	pa = PA(options, PA.AREA + PA.CENTROID, table, MinCellSize, MaxCellSize, MinCirc, MaxCirc)
	#1 -> shows the resulting image | 0-> does not
	pa.setHideOutputImage(1)
	pa.analyze(image)
	
	#get count from table	
	count = table.getCounter()
	
	#catch resulting image
	out = pa.getOutputImage()
	ip = ImageConverter(out)
	ip.convertToRGB()
	# get red slice (to merge with original image for resultsimage)
	red		= getStacks(out,4)
	redout 	= ImagePlus("redoutline",red)

	#adds count to specific (key=filename) to Dictionary
	CountDictionary[filename].append(count)
	calc = IC()

	#cover preprocessed images to give review options
	invertImage(redout)
	#merge redslice and original image for reviewImage
	reviewImage = calc.run("ADD create RGB", imp, redout)

	#reset table to avoid wrong results (this version is a workaround)
	table.reset()
    
	#save image if resultOutput == True
	saveMergedImage(reviewImage,filename)
	
	#return reviewImage

	
#################################################
#### specialized particle Analyser functions ####
#################################################

#Particle Analyser for training (no final image is created)
def PA_training(image,table):


	options = PA.SHOW_ROI_MASKS \
		+ PA.INCLUDE_HOLES  \
	
	roim = RoiManager(True)
	#added PA.CENTROID to get centroids of all particles (X,Y Values)
	pa = PA(options, PA.AREA + PA.CENTROID, table, MinCellSize, MaxCellSize, MinCirc, MaxCirc)
	#hides result image (1)
	pa.setHideOutputImage(1)
	pa.analyze(image)	
	#get count from table	
	count = table.getCounter()
	##load global X,Y-Values for comparison purpose
	global 	PA_list_X
	global	PA_list_Y
	PA_list_X 	= table.getColumn(6)
	PA_list_Y	= table.getColumn(7)
	
	#exception handling (lists are empty)
	if PA_list_X is None:
		PA_list_X = []
	if PA_list_Y is None:
		PA_list_X = []

	#create Array which stores all minimum distances (if smaller than CellDiameter)
	global PA_list_Distance
	PA_list_Distance 	= [CellDiameter] * count
	for i in xrange(0,len(PA_list_Distance)):
		for j in xrange(i+1,len(PA_list_Distance)):
			PA_list_Distance[i] = min(distanceP1P2(PA_list_X[i],PA_list_Y[i],PA_list_X[j],PA_list_Y[j]),PA_list_Distance[i])
	
	table.reset()

#create mask for all found particles (outlines of particles)
def PA_Mask(image,table,minSize):

	options = PA.SHOW_ROI_MASKS \
		+ PA.INCLUDE_HOLES  \
		+ PA.SHOW_MASKS 
		
	roim = RoiManager(True) 
	pa = PA(options, PA.AREA, table, minSize, Double.POSITIVE_INFINITY)
	#hide output image (1)
	pa.setHideOutputImage(1)
	pa.analyze(image)

	#catch resulting image
	out = pa.getOutputImage()
	table.reset()
	return out

#get large particles
def PA_getLargeCount(image):

	options = PA.SHOW_ROI_MASKS \
			+ PA.INCLUDE_HOLES  \

	global oldAvPartSize
	table = ResultsTable()	
	roim = RoiManager(True) 
	pa = PA(options, PA.AREA, table, 0, Double.POSITIVE_INFINITY)
	#hide output image (1)
	pa.setHideOutputImage(1)
	pa.analyze(image)
	
	#if all particles are desolved we return 0
	count = table.getCounter()
	return count

#get too small particles	
def PA_getSmaller(image):

	options = PA.SHOW_ROI_MASKS \
		+ PA.INCLUDE_HOLES  \
		+ PA.SHOW_MASKS 
		
	table = ResultsTable()	
	roim = RoiManager(True) 
	pa = PA(options, PA.AREA, table, MinCellSize, MaxCellSize)
	#hide output images (1)
	pa.setHideOutputImage(1)
	pa.analyze(image)

	#catch resulting image
	out = pa.getOutputImage()
	table.reset()
	return out

#############################
##### save merged image #####
#############################

def saveMergedImage(image,name):
	if resultOutput == True:
		IJ.save(image, outputFolder + "result_" + name )
	else:
		return


###########################
##### image converter #####
###########################

#source : http://rsb.info.nih.gov/ij/developer/api/ij/process/ImageConverter.html

# converts images to 8 bit Thresholding ...
def convertTo8Bit(image):
	ip = ImageConverter(image)
	ip.convertToGray8()
	return image

############################################
##### get MaxEntropy Calculation Value #####
############################################

###############################
##### Threshold-Function ######
###############################

# !!!! IMPORTANT !!!! => This function just works for images with 8Bit depth, don't expect it to work for other types of images 
#(otherwise alternate it the way you want based on the following source code : https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java)

def prepThreshold(image,method,noBlack=False,noWhite=False):


	ip = image.getProcessor()


	data = ip.getHistogram()

	#in case blacks or whites should be left out (usually not needed)
	if noBlack:
		data[0] = 0
	if noWhite:
		data[len(data)-1]=0

	minbin = -1
	maxbin = -1

	for i in xrange(len(data)):
		if data[i] > 0:
			maxbin = i

	j = len(data)-1

	while j >= 0:
		if data[j] > 0:
			minbin = j
		j -= 1


	datLength = maxbin-minbin+1

	#data2 saves just all values, which are represented in the image (cuts off values which where counted 0 times)
	data2 = [0]*datLength

	k = minbin
	while k <= maxbin:
		data2[k-minbin] = data[k]
		k += 1	

	# apply threshold function
	if len(data2) < 2:
		threshold = 0
	elif method == "MaxEntropy":

		threshold = MaxEntropyValue(data2)
	else:
		print "Wrong input"
	
	threshold += minbin
	
	return threshold

################################
##### MaxEntropy-Function ######
################################

def MaxEntropyValue(data):
	#This function was taken from github (https://github.com/fiji/Auto_Threshold/blob/master/src/main/java/fiji/threshold/Auto_Threshold.java)
	#it is included in the Auto_Threshold.java class plugin and can be used with >>>IJ.setAutoThreshold(image,"MaxEntropy dark")<<< 
	#(slightly more sensitive version which can be obtained with calculating the maxentropy value und adjust it by "x")
	#use IJ.setAutoThreshold(image,>>>Value<<<) to get the results you need
	
	###### local variables (types not needed in python/jython) #######
	# "integer" variables
	threshold 	=   -1
	#ih 			= 	0
	#it 			=	0
	lenPix		=	len(data)
	first_bin 	= 	0
	last_bin 	=	0
	# "double" variables
	tot_ent 	=	0.0 # total entropy
	max_ent 	=	0.0 # max entropy
	ent_back 	=	0.0 # entropy of the background pixels at a given threshold
	ent_obj 	=	0.0 # entropy of the object pixels at a given threshold
	# "double-array" variables
	norm_histo 	= 	[0.0] * lenPix  # normalized histogram
	P1 			=	[0.0] * lenPix	# cumulative normalized Histogram
	P2 			=	[0.0] * lenPix  # 
	### end local variables ###
	
	total 		=	0
	for i in xrange(lenPix):
		total+=data[i]
	
	# Important Note: float values needed => danger of obtaining incorrect results 
	for i in range(lenPix):	
		norm_histo[i] = float(data[i])/float(total)

	P1[0] = norm_histo[0]
	# 1.0 to keep the double format active
	P2[0] = 1.0-P1[0]
	for i in xrange(1,lenPix):
		P1[i] = P1[i-1] + norm_histo[i]
		P2[i] = 1.0 - P1[i]

	for i in xrange(lenPix):
		if not(math.fabs(P1[i]) < 2.220446049250313e-16):
			first_bin = i
			break
			
	#determine the last non-zero bin
	last_bin= lenPix-1
	i2 = lenPix-1

	while i2>=first_bin:
		if not(math.fabs(P2[i2]) < 2.220446049250313e-16):
			last_bin = i2
			break
		i2 -=1

	# smallest value
	max_ent = 2.2250738585072014e-308
	j = first_bin

	while j <= last_bin:
		# entropy of background pixels
		ent_back = 0.0
		for i in xrange(j):
			if data[i] != 0:				
				ent_back -= (norm_histo[i]/P1[j]) * math.log(norm_histo[i]  / P1[j])

		# entropy of object pixels
		ent_obj = 0.0
		
		# assign x the value j+1 (loop counter Variable)
		x = j+1
		
		#equal to: for ( x = j + 1; ih < data.length; ih++ ){
		while x < lenPix:

			if data[x] != 0:
				ent_obj -= (norm_histo[x]/P2[j]) * math.log(norm_histo[x]/P2[j])
	
			#while loop counter	
			x+=1
		### end inner while loop

		#total entropy
		tot_ent = ent_back + ent_obj

		
		if max_ent < tot_ent:
			max_ent = tot_ent
			threshold = j
			
		#outer while loop counter
		j += 1

	### end outer while loop

	#return threshold value
	return threshold


###########################################
##### Use Threshold on Image function #####
###########################################


def iteration_ThresHold_init(image):
	
	#load all global variables
	global orgTH 
	global newTH
	global storeTH
	global distance_TH
	global adjust
	global iteration_step_size
	
	#Convert image to 8-bit to make it suitable for analysis
	convertTo8Bit(image)
	
	#obtain threshold
	ownMaxEntropy = prepThreshold(image,"MaxEntropy")
	
	
	#apply adjust to threshold
	storeTH = float(ownMaxEntropy * adjust)
	#####################################################
	#get distance to 255 (highest possible threshold value) // initializing further calculations
	distance_TH 		= float(255-storeTH)
	iteration_step_size 	= float(float(distance_TH/iteration_steps))
	#print "Size of iteration steps: " + str(iteration_step_size)
	ip = image.getProcessor()
	ip.setThreshold(storeTH,255,ImageProcessor.NO_LUT_UPDATE)
	IJ.run(image,"Convert to Mask","")

	#save orgTH
	orgTH = storeTH

#increment threshold by one step
def iteration_ThresHold_increment(image):	
	global newTH
	global storeTH
	
	convertTo8Bit(image)
	#calculate new threshold (float value)
	storeTH = float(storeTH) + float(iteration_step_size)
	#get new to be applied th (ceiled value of stored
	newTH = math.ceil(storeTH)
	#print "Applied TH: " + str(newTH) + " | Stored TH : " + str(storeTH)
	ip = image.getProcessor()
	ip.setThreshold(newTH,255,ImageProcessor.NO_LUT_UPDATE)
	IJ.run(image, "Convert to Mask","")
#### end Thresholding functions


###############################
#### Cell cluster analyser ####
###############################

	#needed:
	# 	1. copied image before threshold is applied
	#	2. MaxCellSize and MinCellSize
	#   3. ~ Threshold startvalue

	#what is the function suppose to do?
	#	-> preprocess all particles bigger than MaxCellSize 
	#How?
	#	By Increasing the threshold value by x and check if the threshold reduced the average particle size enough to be smaller than the MaxCellSize but bigger than 1/2 MaxCellSize

	#what does the function return?:
	# 	The function returns a binary image with all melted off cell clusters (avoiding to return other particles)

def cellClusterAnalyser(image):
	calc = IC()

	#initialize 
	#(needs to be called once seperately to apply first normal Threshold)
	afterSliceImage = image.duplicate()
		
	#apply maxentropy threshold
	iteration_ThresHold_init(afterSliceImage)

	#preprocess for large particles
	LargeCellMask = PA_Mask(afterSliceImage,StatTable,MaxCellSize)
	#preprocess for small particles
	SmallCellMask = calc.run("Substract create", afterSliceImage, LargeCellMask)
	#count left particles
	largeCellCount = PA_getLargeCount(LargeCellMask)

	#close still open images
	afterSliceImage.close()
	LargeCellMask.close()
	#initialization finished

	#continuing with repeating loop
	while largeCellCount != 0:

		#get duplication of imp again
		afterSliceImage = image.duplicate()
		
		#apply TH to sliced image
		iteration_ThresHold_increment(afterSliceImage)
						
		#create Mask for large particles
		LargeCellMask = calc.run("Substract create",afterSliceImage,SmallCellMask)
		#apply watershed (cut almost sparated particles)
		watershed(LargeCellMask)
		#filter all now smaller particles from image and add them to SmallCellMask
		newSmallCellMask = PA_getSmaller(LargeCellMask)
		# merge SmallCellMask with newSmallCellMask
		SmallCellMask = calc.run("add create", SmallCellMask,newSmallCellMask)
		# reduce LargeCellMask by all known small particles
		LargeCellMask = calc.run("Substract create", LargeCellMask, SmallCellMask)
		largeCellCount = PA_getLargeCount(LargeCellMask)

		#close open images
		newSmallCellMask.close()
		LargeCellMask.close()
	
	return  SmallCellMask

#_____________________________________________________ collective function-call-function __________________________________________________
		
###############################
##### Ref_processing Calc #####
###############################
'''
# main image preprocessing and further anylsis steps
# Was used during early development (processing was separated in the later stages of the development)

def Ref_processing(imp):
		
	#______________________________ preprocessing - beginning ______________________________
	#duplicate image to keep original image file
	filename = imp.title
	image = imp.duplicate()
	
	###### inverting image ######
	invertImage(image)
	
	#image.show()
	
	#red,green 	= getStacks(image,3)
	red,green,blue 	= getStacks(image,1)
	greenSlice 	= ImagePlus("GreenSlice_"+filename,green)
	redSlice 	= ImagePlus("RedSlice_"+filename,red)
	blueSlice 	= ImagePlus("BlueSlice_"+filename,blue)
		
	##### substract blue from green #####
	#create object ImageCalculator object (defined as "as IC" at the top)
	calc = IC()

	resultImage = calc.run("Substract create", greenSlice, redSlice)
	#resultImage = calc.run("Add create", resultImage,blueSlice)
	#resultImage2.show()
	
	################################################################################################################################################################################

	# smooth image to get less noise
	image = ImagePlus(filename,filterImage(resultImage,20))
	
	#print "image type " + filename + " " + str(type(image)) + ", imp type " + str(type(imp))
	#______________________________ preprocessing - end ______________________________
	#start analysing process
	analysis_processing(image,filename,imp)

'''

#######################################
##### prepocessed_processing Calc #####
#######################################

#if already preprocessed images need to be analysed
	
def pre_processed(imp):
	
	#prepare for main pre_processing
	filename = imp.title
	
	image = imp.duplicate()
	
	#print "image type " + filename + " " + str(type(image)) + ", imp type " + str(type(imp))
	#start main pre_processing
	analysis_processing(image,filename,imp)

########################################
##### main_anlysis_processing Calc #####
########################################
# calls analysis functions

def analysis_processing(image,filename,imp):
	#__________________ evaluation / validation functions ______________________________________________
	image = cellClusterAnalyser(image)
	
	#watershed for seperation
	watershed(image)

	#print adjust

	#analyze particles 
	#check first if training is initialized
	if param_search:
		PA_training(image,StatTable)
		analyse_lists()

	else:
		#image - current binary image
		# imp  - original image
		ParAna(image,StatTable,CountDic,filename,imp)
	#close all open images
	image.close()
	imp.close()


#____________________________________________________ PART 2 Manual-Auto-Comparison_________________________________________________________	
######################################
##### Part 2 - Comparison functions###
##### Distance function	##############
######################################

#calc distance between manual particle (X1,Y1) and auto particle (X2,Y2)
def distanceP1P2(X1,Y1,X2,Y2):
	a = X1-X2
	b = Y1-Y2
	c = pow(a,2) + pow(b,2)
	return c
	
#######################################
### 		xml-reader				###
#######################################
# For manual counted images file (stored )
# xml file name: CellCounter_result_ + imagename + .xml	
def xml_reader(xmlPath,imageName):
		imageXMLmatch = False
		#delete ending of imageName (.jpg/.tiff etc) -> split string at "."
		#filename.endswith(".JPG") or filename.endswith(".jpg")
		if imageName.endswith(".JPG"):
			imageName = imageName.split(".JPG")[0]
		elif imageName.endswith(".jpg"):
			imageName = imageName.split(".jpg")[0]
		elif imageName.endswith(".TIFF"):
			imageName = imageName.split(".TIFF")[0]
		elif imageName.endswith(".tiff"):
			imageName = imageName.split(".tiff")[0]
		elif imageName.endswith(".TIF"):
			imageName = imageName.split(".TIF")[0]
		elif imageName.endswith(".tif"):
			imageName = imageName.split(".tif")[0]


		
		#save repetetive checks in case of param_search, by adding values to a xmlDic
		if param_search:
			#save all already tested images in dictionary (xmlDic)
			#testing if imageName already exists in xmlDic
			if imageName in xmlDic:
				#read stored values
				imageXMLmatch 	= xmlDic[imageName][0]
				filename 		= xmlDic[imageName][1]
				if imageXMLmatch:
					currentFile = open(os.path.join(xmlPath,filename),"r")
					#pass file to specific analyse function
					XYxmlParser(currentFile)
			#if no entry for current image yet
			else:
				for filename in os.listdir(xmlPath):
					if filename.endswith(".xml") and imageName in filename:
						xmlDic[imageName] = [True,filename]
						imageXMLmatch = True
						currentFile = open(os.path.join(xmlPath,filename),"r")
						#pass file to specific analyse function
						XYxmlParser(currentFile)

		#Iterate through all xml files and check if one matches the image
		#if so extract values
		else:
			for filename in os.listdir(xmlPath):
				if filename.endswith(".xml") and imageName in filename:
					# set xmlDic[imageName] = True and imageXMLmatch = True
					xmlDic[imageName] = True
					imageXMLmatch = True
					currentFile = open(os.path.join(xmlPath,filename),"r")
					#pass file to specific analyse function
					
					XYxmlParser(currentFile)
		if imageXMLmatch == False:
			print "No matching training file for " + imageName
			if "FIE_xml" not in error_dict['FIE']:
				error_dict['FIE'].append("FIE_xml")
			#stop testing for this image
			return 0
		#continue training for this image if matching XML-file was found
		else:
			return 1

#extracts all XY-values from xml file
# basic parsing algorithm (designed for xml files created by imagej cell counter)				
def XYxmlParser(XMLfile):
		#specialised start/end parsing strings (syntax of manual count xml-files)
		inputStart 	= "<Type>1</Type>"
		inputEnd	= "<Type>2</Type>"
		
		#boolean value to indicate if X/Y should be extracted
		parse 		= 0
		
		#with global variables
		global xml_list_X
		global xml_list_Y
		
		#reset global xml_Lists before using
		xml_list_X = []
		xml_list_Y = []
				
		#check for input start
		for line in XMLfile:
			if inputStart in line:
				parse = 1
			#checking for inputEnd -> if so stop loop iteration
			if inputEnd in line:
				parse = 0
				break
			#if parse == True -> parse for X/Y Values
			#keep counter running to know how many values are stored already
			#   <MarkerX>606</MarkerX>
            #	<MarkerY>1610</MarkerY>
			
			#separate different elements from line
			if parse:
				if "<MarkerX>" in line:
					cString = line.split('>',1)
					cString = cString[1].split('<',1)
					xml_list_X.append(int(cString[0]))
					#store x value in xList
				if "<MarkerY>" in line:
					#store y value in yList
					cString = line.split('>',1)
					cString = cString[1].split('<',1)
					xml_list_Y.append(int(cString[0]))

		#close xml-file
		XMLfile.close()		
		
		
		
#PA_Training is implemented next to all other particle analyzer
#######################################
### 	statistical functions		###
#######################################

#get necessary statistical values with manual/auto comparison
#uses global xml/auto_lists to figure amount of matching particles
def analyse_lists():
	#initialize all global arrays
	global xml_list_X,xml_list_Y,manual_hitlist,PA_list_Distance,PA_list_X,PA_list_Y
	#local necessary values (avoid redundant calculations)
	n = len(PA_list_X)
	m =	len(xml_list_X)
	#Catch exception in case nothing could be found in XML-file
	if m == 0 :
		print "XML-file for " + currentImageName + " is corrupt or settings result in no matches!! "
		print "Check if necessary format is correct (Tip: Use the correct Type in Cell-Counter!)"
		if param_search:
			broken_dict[iteration,currentImageName] = [MinCirc,iteration_steps,adjust,MinCellSize]
	elif n == 0:
		print "Error in Particle Analyser ... "
		print "Check your settings, they might lead to no matching particles!"
		print "Continuing with next image..."
		#if validation -> save the broken configuration and fill all values with zeroes
		if param_search:
			broken_dict[iteration,currentImageName] = [MinCirc,iteration_steps,adjust,MinCellSize]
			iterationDic[iteration,currentImageName] = [MinCirc,iteration_steps,adjust,MinCellSize,0,0,0,0,0,0,0]
	else:
		#initialise manual_hitlist -> array size length(xml_list_X), filled with zeros
		manual_hitlist = [0]*m

		#test for every particle from PA_lists if they have a matching particle in xml_lists
		#outer loop - iteration through automatically counted particles
		for i in xrange(0,n):
			#inner loop - iteration through manually counted particles
			for j in xrange(0,m):
				#if particle not hit yet
				if manual_hitlist[j] == 0:
					#calc distance between manual particle in position j and auto particle in position i
					#distanceP1P2(X1,Y1,X2,Y2):
					if distanceP1P2(PA_list_X[i],PA_list_Y[i],xml_list_X[j],xml_list_Y[j]) <= pow(PA_list_Distance[i],2):
					#if distanceP1P2(PA_list_X[i],PA_list_Y[i],xml_list_X[j],xml_list_Y[j]) <= PA_list_Distance[i]:

						#if true set hitlist to 1 and jump to next iteration for outer loop
						manual_hitlist[j] = 1
						break
		#standard values
		hits = 0
		for i in manual_hitlist:
			if i:
				hits = hits+1
		difference = m-hits
		sensitivity = round(float(float(hits)/float(m)),5)	
		wrongCounts = n-hits
		#avoid zero devision issue, in case there was no hit
		if hits == 0:
			hits = 1
		#count error is represented by manual/autocount
		countError 	= math.fabs(round(1.0-(float(m)/float(n)),3))
		precision 	= round(float(float(hits)/float(hits+wrongCounts)),5)
		
		#if we "train", we need 2 key values to differentiate between the iterations and their settings and the images
		if param_search:
			iterationDic[iteration,currentImageName] = [MinCirc,iteration_steps,adjust,MinCellSize,m,n,hits,difference,wrongCounts,countError,sensitivity,precision]
			
#get best setting
#function returns a list with "best" settings and their results
#out best settings are represented by the f-beta-score

#value for correct positives devided by actual positives
def recallScore(iter,name):
	#recall = TP/(TP+FN) 
	#TP -> hits -> iterationDic[iter,6]
	#FN -> missed -> iterationDic[iter,7]
	recall = float(iterationDic[iter,name][6])/float((float(iterationDic[iter,name][6])+float(iterationDic[iter,name][7])))	
	return recall 

#value for correct positives devided by predicted positives
def precScore(iter,name):
	#precision = TP/(TP+FP) 
	#TP -> hits -> iterationDic[iter,6]
	#FP -> wrong counts -> iterationDic[iter,8]
	precision = float(iterationDic[iter,name][6])/float((float(iterationDic[iter,name][6])+float(iterationDic[iter,name][8])))
	return precision
	
	
#calculates the F-Score for all images in an iteration
#The F-Score represents a measurement of accuracy (harmonic mean of precision and recall)
def f1Score():
	qualityDic.clear()
	qCountDic.clear()
	
	for iter,name in iterationDic:
		#add the calculated value to the specific qualityDic items (based on the iteration)
		#f1 score = 2* ((Precision*Recall)/(Precision+Recall))
		#if there are no hits, there will be no precision and no recall value (they will be equal to 0) -> f1 score is 0
		if (iterationDic[iter,name][8] == 0):
			f1			= 0
		else:
			#get precision and recall for specific image and iteration
			precision 	= precScore(iter,name)
			recall 		= recallScore(iter,name)
			f1			= (1+float(pow(f_beta,2)))* float(float(precision * recall)/float(float(pow(f_beta,2))*precision + recall))
		#f1			= 2* float(float(precision * recall)/float(precision + recall))	
		#if dictionary has no entries yet:
		#handles exception if no values have been saved yet for our iteration key
		qualityDic.setdefault(iter,(0))
		qCountDic.setdefault(iter,(0))
		#adds values to the dictionary
		qualityDic[iter] = qualityDic[iter] + f1
		qCountDic[iter] += 1
		#add F1 Score to image in iteration dictionary
		iterationDic[iter,name].append(f1)
	#now devide score saved in qualityDic[iter][0] by the count of fitting images
	for iter in qualityDic:
		qualityDic[iter] = float(qualityDic[iter])/float(qCountDic[iter])		
	#return values are the maximum and minimum value in our qualityDic
	#return max(qualityDic, key=qualityDic.get),min(qualityDic,key=qualityDic.get)


#Filter function for images - marks images as not usable if wanted by user
#The value, which determines if a image is "broken", is defined by "accuracyTH"
def brokenImageFilter():
	
	for iter,name in iterationDic:
		if (iterationDic[iter,name][6] == 0):
			f1			= 0
		else:
			#get precision and recall for specific image and iteration
			precision 	= precScore(iter,name)
			recall 		= recallScore(iter,name)
			f1			= 2* float(float(precision * recall)/float(precision + recall))	
		#if the countError (index 5) exceeds our accuracyTH, copy it to broken_dict and delete it from iterationDic
		if f1 < accuracyTH:
			#create entry for broken_dict
			#copy over values (minCirc,iteration_steps,adjust)
			broken_dict[iter,name] = [iterationDic[iter,name][0],iterationDic[iter,name][1],iterationDic[iter,name][2],iterationDic[iter,name][3],	f1]
	#delete all entries with a lower score than TH from iterationDic and add the names to brokenIterNameDic
	for iter,name in broken_dict:
		del iterationDic[iter,name]
		brokenIterNameDic.setdefault(iter,[])
		brokenIterNameDic[iter].append(name) 

#calculates maxCellSize and CellDiameter from MinCellSizeDialogListener
def CellDia(maxC):

	#circle area  = Pi * r^2, we increase it by 5% (most cells are not perfectly round -> we get more accurate hit results)
	#that means we can get the radius from maxC: sqrt(maxC/pi) -> diameter 2*r
	pi 		= math.pi
	diaC 	= 2*sqrt(maxC/pi) * 1.05 #increased by 5% and round to 2 float points
	
	return diaC
	
	
#______________________________________________ Function-Area-end_________________________________________


#####################################
##### Interface calls - Dialogs #####
#####################################

# Dialog Options : http://rsb.info.nih.gov/ij/developer/api/ij/gui/GenericDialog.html

# used to get folder for all xml files (manual count XY values)
def callXMLDialog():

	#____choose origin folder for files ___
	dialog = GenericDialog("Cell Image Analyser")
	dialog.addMessage("Please choose the path to the XML-folder-location.")
	dialog.addMessage("")
	dialog.showDialog()  
	
	if dialog.wasCanceled():
		print "User canceled dialog!"  
		return  
	
	return True
	
#XML dialog
def xml_folder_chooser():
	
	direc = DirectoryChooser("Choose folder with manual counted xml-files.")
	XMLfolder = direc.getDirectory()	
	
	if folder is None:
		print "User canceled the dialog! Exiting ..."
	else:
		print "Selected Folder :", XMLfolder	

		return XMLfolder	


def callImageDialog():

	#____choose origin folder for files ___
	dialog = GenericDialog("Cell Image Analyser")
	dialog.addMessage("           Welcome to the Cell Image Analyser.")
	dialog.addMessage("")
	dialog.addMessage("Please start by choosing the folder containing your images.")
	dialog.addMessage("")
	dialog.showDialog()  
	
	if dialog.wasCanceled():
		print "User canceled dialog!"  
		return  
	
	return True



def image_folder_chooser():

	#____chose input/output_Folder ___

	direc = DirectoryChooser("Choose folder with images.")
	folder = direc.getDirectory()

	if folder is None:
		print "User canceled the dialog! Exiting ..."
	else:
		print "Selected Folder :", folder	
		
		#create dialog object
		dialog = GenericDialog("Cell Image Analyser")  

		dialog.addMessage("Analysing Options:")
		types = [ "Automatic cell counting" , "Best Parameter Search"]
		#default choice  = basic
		dialog.addChoice("Choose:", types, types[0])  
		dialog.showDialog() 
		
		if dialog.wasCanceled():
			print "User canceled dialog!"  
			return  

		settingChoice 	= dialog.getNextChoice()

		return settingChoice, folder	
		
		
################################
##### call advanced Dialog #####
################################


def advancedDialog():

		#create dialog object
		dialog = GenericDialog("Cell Image Analyser - evaluation")  
		#create input options  
		dialog.addStringField("Analysed by:", "annonymous")
		dialog.addStringField("Output postfix:","standard") 
		dialog.addCheckbox("Save review images?", True)
		dialog.addNumericField("Minimum cell size (Pixel):", 30 , 1)
		dialog.addNumericField("Maximum cell size (Pixel):", 150 , 1)
		dialog.addNumericField("Minimum circularity:", 0.3, 1)
		dialog.addNumericField("Maximum circularity:", 1.0, 1)
		dialog.addNumericField("Number of iterations to saturate threshold:",20,0)
		dialog.addNumericField("Adjust Threshold in (+-)%:", 5,1)

		dialog.showDialog()  
		#  
		if dialog.wasCanceled():
			print "User canceled dialog!"  
			return  
	  	# Read out the values 
	 	creator			= dialog.getNextString()
		outputPostfix	= dialog.getNextString()
		
		resultOutput	= dialog.getNextBoolean()	
		MinCellSize 	= dialog.getNextNumber()
		MaxCellSize 	= dialog.getNextNumber()
		MinCirc			= dialog.getNextNumber()
		MaxCirc			= dialog.getNextNumber()
		iteration_steps			= dialog.getNextNumber()
		adjust 			= dialog.getNextNumber()
		
		
		#MaxCellSize and CellDiameter are based on MinCellSize
		CellDiameter = CellDia(MaxCellSize)
		#exception handling: Min / MaxCirc
		#In case circularity definition is unclear
		#return creator, outputPostfix, resultOutput, fileInput, MinCellSize, MaxCellSize, MinCirc, MaxCirc, iteration_steps, adjust, CellDiameter
		return creator, outputPostfix, resultOutput, MinCellSize, MaxCellSize, MinCirc, MaxCirc, iteration_steps, adjust, CellDiameter

		return
		


#changed training dialog (MinCellSize instead Contrast)
def param_search_dialog():
		#create dialog object
		#This dialog is used to adjust training range for usefull values
		dialog = GenericDialog("Cell Image Analyser - statistics")  
		
		#create input options
		dialog.addStringField("Analysed by:", "annonymous")
		dialog.addStringField("Result prefix:","")

		#Static values
		dialog.addMessage("Static values:")
		dialog.addNumericField("F-score beta value:", 1, 1)
		dialog.addNumericField("MaxCellSize (Pixel):",300,0)
		dialog.addNumericField("MaxCirc:", 1.0, 2)
		
		#outlier filter
		dialog.addCheckbox("Apply accuracy filter (filter by F-Score)?",False)
		dialog.addNumericField("Exclude image, if F-Score < ", 0.85, 2)
		#preference for suggested settings
		pref = ["Maximum number of images", "Highest accuracy (F1-Score)"] 
		dialog.addChoice("Sort results by:", pref, pref[0])

		dialog.addMessage("Dynamic values:")
		
		#dynamic values
		dialog.addCheckbox("Keep minimum circularity static?",False)
		dialog.addNumericField("Lowest value minimal circularity:", 0.3,2)
		dialog.addNumericField("Incrementation step size minimum circularity:", 0.2,2)
		dialog.addNumericField("Maximum value for minimum circularity:", 0.8,2)


		dialog.addCheckbox("Keep number of iterations static?",False)
		dialog.addNumericField("Lowest number iterations:", 10,0)
		dialog.addNumericField("Incrementation step size iterations:", 5,0)
		dialog.addNumericField("Maximum number iterations:", 30,0)
		
		dialog.addCheckbox("Keep threshold adjustment static?",False)
		dialog.addNumericField("Start value threshold adjustment (+-%):", 0,0)
		dialog.addNumericField("Incrementation step size threshold adjustment (+%):", 5,0)
		dialog.addNumericField("Maximum threshold adjustment (+-%):", 15,0)

		#MinCellSize as a not static parameter (MaxCellSize and Diameter are directly linked to MinCellSize)
		dialog.addCheckbox("Keep minimum cell size static?",False)
		dialog.addNumericField("Lowest minimum cell size (pixel):", 30, 0)
		dialog.addNumericField("Incrementation step size minimum cell size?:",  50,0)
		dialog.addNumericField("Largest minimum cell size (pixel):",150,0)




		dialog.showDialog()
		
		#  
		if dialog.wasCanceled():
			print "User canceled dialog!"  
			return  
	  	# Read out the values 
		creator				= dialog.getNextString()
		outputPrefix		= dialog.getNextString()

		# fscore value
		f_beta				= dialog.getNextNumber()
		MaxCellSize			= dialog.getNextNumber()
		MaxCirc				= dialog.getNextNumber()
		
		#outlier choice
		outlierChoice		= dialog.getNextBoolean()
		accuracyTH			= dialog.getNextNumber()

		#Preference for best setting
		prefSetting			= dialog.getNextChoice()

		#dynamic values 
		MinCirc				= dialog.getNextNumber()
		MinCircStep			= dialog.getNextNumber()
		MaxMinCirc			= dialog.getNextNumber()

		minStatic			= dialog.getNextBoolean()
		iteration_steps		= dialog.getNextNumber()
		iteration_inc		= dialog.getNextNumber()
		max_iterations		= dialog.getNextNumber()
		iteration_static	= dialog.getNextBoolean()
		adjust 				= dialog.getNextNumber()
		adjustStep			= dialog.getNextNumber()
		adjustMax			= dialog.getNextNumber()
		adjustStatic		= dialog.getNextBoolean()
		MinCellSize			= dialog.getNextNumber()
		cellStep			= dialog.getNextNumber()
		cellRange			= dialog.getNextNumber()

		
		cellStatic			= dialog.getNextBoolean()


		#Background subtraction


		#MaxCellSize and CellDiameter are based on MinCellSize
		CellDiameter = CellDia(MaxCellSize)
		#exception handling: Min / MaxCirc
		#In case circularity definition is unclear
			
		return creator, outputPrefix, MinCirc, MaxMinCirc, MinCircStep, iteration_steps, iteration_inc, adjust, adjustStep, cellStep, prefSetting, MinCellSize, MaxCellSize, MaxCirc, CellDiameter, minStatic, iteration_static, adjustStatic, cellStatic, outlierChoice, accuracyTH, cellRange, max_iterations,  adjustMax, f_beta

														##################################
#_______________________________________________________##### File-writing functions #####_______________________________________________________
														##################################


######################################
##### write counting-result file #####
######################################

def saveCountResults(CountDic):
	#open txt file
	if resultOutput==False:
		print "Saving results in :", folder
		results = open(folder + "Evaluation_Result.txt","w")

		#write timestamp and creator
		ctime = datetime.now()
		results.write("Analysis done by	" + creator + "\n")
		results.write("Date and time	" + str(ctime.day) + "." + str(ctime.month) + "." + str(ctime.year) + "\t" + str(ctime.hour) + ":" + str(ctime.minute) + ":" + str(ctime.second) + "\n")
		results.write("Pub version		" + str(pubVersion) + "\n")
		results.write("Feature version	" + str(Feature_Version) + "\n")
		results.write("Cell Counter version	" + str(CellCounterVersion) + "\n")
		results.write("Modus	" + mode + "\n")
		#occured errors
		results.write("____________________________________" + "\n")
		results.write("List of occured errors:\n")
		results.write("\n")
		value_counter = 0
		for item in error_dict:
			value_counter = value_counter + len(error_dict[item])
			
		if value_counter == 0:
			results.write("No errors detected.\n")
		else:
			for item in error_dict:
				if len(error_dict[item]) != 0:
					results.write(item + ":	")
					for val in error_dict[item]:
						results.write(val + "	")
					results.write("\n")
		#write used values to file
		results.write("____________________________________" + "\n")
		results.write("MinCellSize	" + str(MinCellSize) + "\n")
		results.write("MaxCellSize	" + str(MaxCellSize) + "\n")
		results.write("MinCirc	" + str(MinCirc) + "\n")
		results.write("MaxCirc	" + str(MaxCirc) + "\n")
		results.write("CellDiameter	" + str(CellDiameter) + "\n")
		results.write("Threshold iterations	" + str(iteration_steps) + "\n")
		results.write("Threshold adjusted	" + str(adjust) + "\n")

 
		results.write("____________________________________" + "\n")
		results.write("Images tested	" + str(globalImageCount) + "\n")
		results.write("calculation-time	" + str(finalTime) + "\n")
		results.write("Average time per image	" + str(finalTime/globalImageCount) + "\n")
		results.write("____________________________________" + "\n")
			
		results.write("\n")
		results.write("Image-name \t count \n")
		results.write("\n")
		#writes all items to file
		for item in CountDic:
			
			results.write(repr(item).replace('u','') + "\t")
			clist = CountDic[item]
			for i in range(len(clist)):
				results.write(str(clist[i])+"\t")
			results.write("\n")
		results.close()
		
	else: 
		print "Saving images and results in :", outputFolder
		results = open(outputFolder + "Evaluation_Result.txt","w")


		#write timestamp
		ctime = datetime.now()
		results.write("Analysis done by	" + creator + "\n")
		results.write("Date and time	" + str(ctime.day) + "." + str(ctime.month) + "." + str(ctime.year) + "\t" + str(ctime.hour) + ":" + str(ctime.minute) + ":" + str(ctime.second) + "\n")
		results.write("Pub version		" + str(pubVersion) + "\n")
		results.write("Feature version	" + str(Feature_Version) + "\n")
		results.write("Cell Counter version	" + str(CellCounterVersion) + "\n")
		results.write("Modus	" + mode + "\n")
		#occured errors
		results.write("____________________________________" + "\n")
		results.write("List of occured errors:\n")
		results.write("\n")
		value_counter = 0
		for item in error_dict:
			value_counter = value_counter + len(error_dict[item])
			
		if value_counter == 0:
			results.write("No errors detected.\n")
		else:
			for item in error_dict:
				if len(error_dict[item]) != 0:
					results.write(item + ":	")
					for val in error_dict[item]:
						results.write(val + "	")
					results.write("\n")

		#write used values to file
		results.write("____________________________________" + "\n")
		results.write("MinCellSize	" + str(MinCellSize) + "\n")
		results.write("MaxCellSize	" + str(MaxCellSize) + "\n")
		results.write("MinCirc	" + str(MinCirc) + "\n")
		results.write("MaxCirc	" + str(MaxCirc) + "\n")
		results.write("CellDiameter	" + str(CellDiameter) + "\n")
		results.write("Threshold iterations	" + str(iteration_steps) + "\n")
		results.write("Threshold adjusted	" + str(adjust) + "\n")

 
		results.write("____________________________________" + "\n")
		results.write("Images tested	" + str(globalImageCount) + "\n")
		results.write("calculation-time	" + str(finalTime) + "\n")
		results.write("Average time per image	" + str(finalTime/globalImageCount) + "\n")
		results.write("____________________________________" + "\n")		
		
		results.write("\n")
		results.write("Image-name \t count \n")
		results.write("\n")
		#writes all items to file
		for item in CountDic:
			
			results.write(repr(item).replace('u','') + "\t")
			clist = CountDic[item]
			for i in range(len(clist)):
				results.write(str(clist[i])+"\t")
			results.write("\n")
		results.close()


#iteration file writer 
def saveTrainingResults():
	
	print "Saving results in :", folder
	if outputPrefix == "":
		results = open(folder + "_Optimisation_Results.txt","w")
	else:
		results = open(folder + outputPrefix + "_Optimisation_Results.txt","w")
		
		
	#write timestamp and creator
	ctime = datetime.now()
	results.write("Analysis done by 	" + creator + "\n")
	results.write("Date and time  	" + str(ctime.day) + "." + str(ctime.month) + "." + str(ctime.year) + "\t" + str(ctime.hour) + ":" + str(ctime.minute) + ":" + str(ctime.second) + "\n")
	results.write("Pub version		" + str(pubVersion) + "\n")
	results.write("Feature version	" + str(Feature_Version) + "\n")
	results.write("Cell Counter version	" + str(CellCounterVersion) + "\n")
	results.write("Modus	" + mode + "\n")
	#occured errors
	results.write("____________________________________" + "\n")
	results.write("List of occured errors:\n")
	results.write("\n")
	value_counter = 0
	for item in error_dict:
		value_counter = value_counter + len(error_dict[item])
		
	if value_counter == 0:
		results.write("No errors detected.\n")
	else:
		for item in error_dict:
			if len(error_dict[item]) != 0:
				results.write(item + ":	")
				for val in error_dict[item]:
					results.write(val + "	")
				results.write("\n")
	#write used values to file
	results.write("____________________________________" + "\n")

	results.write("Applied F-Score	" + str(f_beta) + "\n")
	results.write("______\n")

	if conStatic:
		results.write("CellSize was kept static at start value	" + str(start_value_array[3]) + "\n")
		results.write("MaxCellSize	" + str(MaxCellSize) + "\n")
		results.write("CellDiameter	" + str(start_value_array[4]) + "\n")
	else:
		results.write("MinCellSize(start)	" + str(start_value_array[3]) + "\n")
		results.write("Max-MinCellSize	" + str(cellRange) + "\n")
		results.write("CellSize-step	" + str(cellStep) + "\n")
		results.write("MaxCellSize	" + str(MaxCellSize) + "\n")
		results.write("CellDiameter	" + str(start_value_array[4]) + "\n")
	results.write("______\n")	
	if minStatic:
		results.write("Minimum circularity was kept static at start value	" + str(start_value_array[0]) + "\n")
	else:
		results.write("MinCirc(start)	" + str(start_value_array[0]) + "\n")
		results.write("Circularity-step	" + str(MinCircStep) + "\n")
		results.write("MaxMinCirc	" + str(MaxMinCirc) + "\n")	
	results.write("MaxCirc	" + str(MaxCirc) + "\n")
	
	results.write("______\n")
	if iteration_static:
		results.write("Iterations were kept static at start value	" + str(start_value_array[1]) + "\n")
	else:
		results.write("Threshold iterations(start)	" + str(start_value_array[1]) + "\n")
		results.write("Threshold iteration increase(start)	" + str(iteration_inc) + "\n")
		results.write("Maximum number of threshold iterations	" + str(max_iterations) + "\n")
		
	results.write("______\n")	
	if adjustStatic:
		results.write("Adjuster was kept static at start value	" + str(start_value_array[2]) + "\n")
	else:
		results.write("Threshold adjusted(start)	" + str(start_value_array[2]) + "\n")
		results.write("Threshold step-size	" + str(adjustStep) + "\n")
		results.write("Threshold maximum	" + str(adjustMax) + "\n")
	results.write("______\n")

	results.write("____________________________________" + "\n")
	results.write("Iterations	" + str(iteration-1) + "\n")
	results.write("Images tested	" + str(globalImageCount) + "\n")
	results.write("calculation-time	" + str(finalTime) + "\n")
	results.write("Average time per iteration	" + str(finalTime/iteration) + "\n")
	results.write("Average time per Image	" + str(finalTime/globalImageCount) + "\n")
	results.write("Average time image:iteration	" + str(finalTime/globalImageCount/(iteration-1)) + "\n")
	results.write("____________________________________" + "\n")

	results.write("\n")

	#save qualityDic values with key and qCountDic values in list
	#create empty list
	printDic = []
	#put all elements from both dictionaries in list (based on iteration)
	for item in qualityDic:
		printDic.append({'iteration': item ,'F1Score':qualityDic[item],'ImageCount':qCountDic[item]})
		
	bestSetting 	= 0
	worstSetting 	= 0
	repeaterC 		= 0
	
	
	
	if outlierChoice:
		results.write("Iteration	" + "Average F_"+ str(f_beta) +"-Score	" + "Viable Images	"  + "Min Circularity	" + "Number of iterations	" + "Threshold adjustment	" + " Min CellSize	" + "Filtered Images" +  "\n \n")
	else: 
		results.write("Iteration	" + "Average F_"+ str(f_beta) +"-Score	" + "Viable Images	" + "Min Circularity	" + "Number of iterations	" + "Threshold adjustment	" + " Min CellSize	"+ "\n \n")
	if (prefSetting == "Maximum number of images"):
		for key in sorted(printDic, key = itemgetter('ImageCount','F1Score'),reverse=True):
			results.write(str(key['iteration']) +  "	" + str(key['F1Score']) + "	" + str(key['ImageCount']) + "	" + str(settingDic[key['iteration']][0]) + "	" + str(settingDic[key['iteration']][1]) + "	" + str((settingDic[key['iteration']][2]-1)*100) + "	" + str(settingDic[key['iteration']][3]) + "	")
			#get best and worst settings
			if repeaterC == 0:
				bestSetting = key['iteration']
			if repeaterC == len(printDic)-1:
				worstSetting = key['iteration']
			
			if outlierChoice:
				results.write("[")
				for i in xrange(0,len(brokenIterNameDic[key['iteration']])):
					results.write(str(brokenIterNameDic[key['iteration']][i]))
					if i < len(brokenIterNameDic[key['iteration']])-1:
						results.write(", ")
				results.write("]"+"\n")
			else:
				results.write("\n")
				
			repeaterC += 1
		results.write("\n")
		results.write("____________________________________" + "\n")
		results.write("\n")
	else:
		for key in sorted(printDic, key = itemgetter('F1Score','ImageCount'),reverse=True):
			print key
			
			results.write(str(key['iteration']) +  "	" + str(key['F1Score']) + "	" + str(key['ImageCount']) + "	" + str(settingDic[key['iteration']][0]) + "	" + str(settingDic[key['iteration']][1]) + "	" + str((settingDic[key['iteration']][2]-1)*100) + "	" + str(settingDic[key['iteration']][3]) + "	")
			
						
			if repeaterC == 0:
				bestSetting = key['iteration']
			if repeaterC == len(printDic)-1:				
				worstSetting = key['iteration']			
				
			if outlierChoice:
				results.write("[")
				for i in xrange(0,len(brokenIterNameDic[key['iteration']])):
					results.write(str(brokenIterNameDic[key['iteration']][i]))
					if i < len(brokenIterNameDic[key['iteration']])-1:
						results.write(", ")
				results.write("]"+"\n")
			else:
				results.write("\n")
			repeaterC += 1
		results.write("\n")
		results.write("____________________________________" + "\n")
		results.write("\n")
	if outlierChoice:
		results.write("Best found setting  :	" + "Iteration : " + str(bestSetting) + "	F_"+ str(f_beta) +"-Score : " + str(qualityDic[bestSetting])   + "	Valid Images : " + str(qCountDic[bestSetting]) + "	Min Circularity : " + str(settingDic[bestSetting][0]) + "	Number of iterations : " + str(settingDic[bestSetting][1]) + "	Threshold adjustment : " + str((settingDic[bestSetting][2]-1)*100) + "	Min CellSize : " + str(settingDic[bestSetting][3]) + "\n")
		results.write("Worst found setting :	" + "Iteration : " + str(worstSetting) + "	F_"+ str(f_beta) +"-Score : " + str(qualityDic[worstSetting])  + "	Valid Images : " + str(qCountDic[worstSetting]) +"	Min Circularity : " + str(settingDic[worstSetting][0]) + "	Number of iterations : " + str(settingDic[worstSetting][1]) + "	Threshold adjustment : " + str((settingDic[worstSetting][2]-1)*100) + "	Min CellSize : " + str(settingDic[worstSetting][3]) + "\n")
	else:
		results.write("Best found setting  :	" + "Iteration : " + str(bestSetting) + "	F_"+ str(f_beta) +"-Score : " + str(qualityDic[bestSetting])   +  "	Min Circularity : " + str(settingDic[bestSetting][0]) + "	Number of iterations : " + str(settingDic[bestSetting][1]) + "	Threshold adjustment : " + str((settingDic[bestSetting][2]-1)*100) + "	Min CellSize : " + str(settingDic[bestSetting][3]) + "\n")
		results.write("Worst found setting :	" + "Iteration : " + str(worstSetting) + "	F_"+ str(f_beta) +"-Score : " + str(qualityDic[worstSetting])  +"	Min Circularity : " + str(settingDic[worstSetting][0]) + "	Number of iterations : " + str(settingDic[worstSetting][1]) + "	Threshold adjustment : " + str((settingDic[worstSetting][2]-1)*100) + "	Min CellSize : " + str(settingDic[worstSetting][3]) + "\n")
	results.write("____________________________________" + "\n \n")
	results.write("\n")
	results.write("Overview for all calculated values: \n")
	
	results.write("\n")
	results.write("Iteration	" + "Name	" + "MinCirc	" + "iteration_steps	" + "adjust	" + "MinCellSize	"  + "manual	" + "auto	" + "hits	"+ "difference	" + "wrongCounts	" + "countError	" + "sensitivity	"  + "precision	" + "f_"+ str(f_beta) +"-score	" + "\n")
	#writes all items to file
	results.write("\n")
	for iter,name in sorted(iterationDic):
		results.write(str(iter) + "	" + str(name) + "	")
		cList = iterationDic[iter,name]
		#you don't want to know why - monstrosity!
		adjust_find = 0
		for i in cList:
			adjust_find += 1
			if adjust_find == 3:
				results.write(str((i-1)*100) + "\t")
			else:
				results.write(str(i) + "\t")

		results.write("\n")
	
	#write the images and their iteration, we couldnt get results for or which exceeded our maximum values, into file
	#for proper sorting (more dimensional sorting -> https://www.youtube.com/watch?v=avthXghx0a0 )
	#use prefSetting = ["Maximum number of images", "Highest accuracy"] 
	
	results.write("____________________________________" + "\n")
	results.write("broken images and their settings: \n")
	if outlierChoice:
		results.write("\n")
		results.write("Outlier-Filter active, set to " + str(accuracyTH))
	results.write("\n")
	if len(broken_dict) == 0:
		results.write("No broken images/Settings found! \n")
		results.write("\n")
	else:
		for iter,name in sorted(broken_dict):
			results.write(str(iter) + "	" + str(name) + "	")
			cList = broken_dict[iter,name]
			adjust_find = 0
			for i in cList:
				adjust_find += 1
				if adjust_find == 3:
					results.write(str((i-1)*100) + "\t")
				else:
					results.write(str(i) + "\t")
		

	results.close()

	

#########################
####### main-body #######
#########################

print "Starting image-analysis..."


##################################
#______ 	open Dialog and value inputs 	______#

any_input = callImageDialog()
if any_input is not None:
	settingChoice = image_folder_chooser()
	#which mode we used
	mode = ""
	globalImageCount = 0

	if settingChoice is not None:  
		setting, folder = settingChoice  #, helpOutput 


		if setting == "Automatic cell counting":
			dialogInput = advancedDialog()
			advancedChoice 	= True
			param_search 	= False
		else:
			dialogInput = param_search_dialog()
			advancedChoice 	= False
			param_search 	= True
			
			
		if dialogInput is not None:  
			

			#################################
			###### Advanced Settings ########
			#################################
			#returns several values
			if advancedChoice: 
				creator, outputPostfix, resultOutput, MinCellSize, MaxCellSize, MinCirc, MaxCirc, iteration_steps, adjust, CellDiameter = dialogInput
				
				mode = "Advanced Settings"



			
				print "Choice:"
				if resultOutput:
					print "Review images will be created"
				else:
					print "No review images will be created"
				#Advanced Dialog input EXCEPTIONS

				if MinCellSize > MaxCellSize:
					MinCellSize = MaxCellSize
					print "MinCellSize bigger than MaxCellSize! Set to MaxCellSize."
					error_dict['PIE'].append("PIE_CellSize")
				if MinCirc >= MaxCirc:
					MinCirc = MaxCirc - 0.1
					print "MinCirc bigger than MaxCirc! Set to MaxCirc - 0.1."
					error_dict['PIE'].append("PIE_Circ")
				if iteration_steps <= 0:
					iteration_steps = 10
					print "Scaling Threshold value too low! Set to default value 1.05."
					error_dict['PIE'].append("PIE_TH")
				if abs(adjust) > 30:
					print "Threshold-adjustment potentially unreasonable (more than 30% change)!"
					if "PIE_TH" not in error_dict['PIE']:
						error_dict['PIE'].append("PIE_TH")

				print "MinCellSize (Pixel):", MinCellSize
				print "MaxCellSize (Pixel):", MaxCellSize			
				print "MinCirc (Pixel):", MinCirc
				print "MaxCirc (Pixel):", MaxCirc
				print "CellDiameter  :", CellDiameter
				print "Threshold Iterations :", iteration_steps
				print "Calculated Threshold will be adjusted by ", adjust, "%"

				adjust = float(float(1)+ adjust*0.01)
			

				
			###################################
			###### Statistical Algorithm ######
			###################################

			elif param_search:
				creator, outputPrefix, MinCirc, MaxMinCirc, MinCircStep, iteration_steps, iteration_inc, adjust, adjustStep, cellStep, prefSetting, MinCellSize, MaxCellSize, MaxCirc, CellDiameter, minStatic, iteration_static, adjustStatic, conStatic, outlierChoice, accuracyTH, cellRange, max_iterations, adjustMax, f_beta = dialogInput
				#add start-values to list
				start_value_array.extend([MinCirc,iteration_steps,adjust,MinCellSize,CellDiameter])
				
				if outlierChoice:
					mode = "Parameter search with Accuracy-filter"

				else:
					mode = "Parameter search without Accuracy-filter"
				#cover not covert values
				resultOutput = False

				#create necessary dictionaries
				#statistical dictionary (reset after each iteration)
				statDic 		= {}
				#iteration Dictionary saves for every iteration all the image results
				iterationDic 	= {}
				#setting dictionary saves settings for each iteration
				settingDic 		= {}
				#quality Dictionary saves the qualityScore for each iteration
				qualityDic 		= {}
				#qCountDic -> counts images which passed the filter
				qCountDic 		= {}			
				
				#create XML-Dictionary (used to avoid redundant check for XML files)
				xmlDic	= {}
				
				#get location of training files (xml format)
				
				any_input = callXMLDialog()	
				if any_input is not None:
					XMLfolder = xml_folder_chooser()

					print "Training statistical values"
					if outlierChoice:
						if accuracyTH > 1 or accuracyTH == 0:
							print "Accuracy-threshold higher than 100% or 0%! Set to default value 0.85 (85%)."
							accuracyTH = 0.85
							error_dict['PIE'].append("PIE_accuracyTH")

						print "Accuracy-Filter on, set to " + str(accuracyTH) + "."
					else:			
						print "Accuracy-Filter off"

					if f_beta < 0:
						print "F-score beta value < 0! Set to default value 1 (balanced F-score)."
						f_beta = 1
						error_dict['PIE'].append("PIE_F_score")
					print "F-Score beta value	 :", f_beta
					print "Circularity static?      :", minStatic
					print "MinCirc (Pixel)          :", MinCirc
					print "MinCirc steps            :", MinCircStep
					print "MaxMinCirc               :", MaxMinCirc
					print "MaxCirc (Pixel)          :", MaxCirc
					print "Scaler static?           :", iteration_static
					print "Threshold iterations        :", iteration_steps
					print "Iteration increase per step  :", iteration_inc
					print "Maximum number of iteartions:", max_iterations
					print "Adjustment static?       :", adjustStatic
					print "Threshold adjustment     :", adjust
					print "Threshold adjust steps   :", adjustStep
					print "Threshold max adjustment :", adjustMax
					print "MinCellSize static?      :", conStatic
					print "MinCellSize (Pixel)      :", MinCellSize
					print "Max-MinCellSize (Pixel)  :", cellRange
					print "MinCell - steps          :", cellStep
					print "MaxCellSize (Pixel)      :", MaxCellSize	
					print "CellDiameter             :", CellDiameter					
				#create output-folder, if wanted		
			if resultOutput:
				if not outputPostfix:
					outputFolder = folder + "output" + "/"
				else:
					outputFolder = folder + "output " + outputPostfix + "/"
				if not os.path.exists(outputFolder):
					os.makedirs(outputFolder)
			
			#_____			 end dialog					______#
				
			#start calc time (for Benchmarking)
			print "Start Processing ..."
			a = datetime.now()
			
			#___________________________________ "Training" ___________________________________________#
			if param_search and XMLfolder is not None:
				
				#Loop iterations:
				#Here we cover all necessary steps to cover non static parameters
				#step count for all parameters are saved in dynamic array
				stepLst = []

				####################### INPUT EXCEPTIONS ##########################
				#add all values
				#Circularity steps
				if MaxMinCirc > MaxCirc:
					print "MaxMinCirc bigger than MaxCirc -> set to MaxCirc - 0.2! (" + str(MaxCirc-0.2) + ")"
					MaxMinCirc = MaxCirc - 0.
					error_dict['PIE'].append("PIE_Circ") 
				if MinCirc > MaxMinCirc:
					print "MinCirc bigger than MaxMinCirc -> set to MaxMinCirc! (" + str(MaxMinCirc) + ")"
					MinCirc = MaxMinCirc
					if "PIE_Circ" not in error_dict:
						error_dict['PIE'].append("PIE_Circ")
				if MinCirc > MaxCirc:
					print "MinCirc bigger than MaxCirc -> set to MaxMinCirc! (" + str(MaxMinCirc) + ")"
					MinCirc = MaxMinCirc
					if "PIE_Circ" not in error_dict:
						error_dict['PIE'].append("PIE_Circ")
				#lifehack ( +0.01, otherwise sometimes wrong results? e.g. 3.0 -> 2.0)
				stepLst.append(int(floor((MaxMinCirc-MinCirc)/MinCircStep + 0.01)))
				#number of iterations
				if iteration_steps <= 0:
					iteration_steps = 10
					print "Number of iterations too low! Set to default value 10."
					error_dict['PIE'].append("PIE_TH")
				#lifehack ( +0.01, otherwise sometimes wrong results? e.g. 3.0 -> 2.0)
				stepLst.append(int(float(float(max_iterations-iteration_steps)/iteration_inc)+0.01))
								
				#adjust steps
				if abs(adjust) > 80:
					print "Threshold-adjustment potentially unreasonable (more than +- 80% change)!"
					if "PIE_TH" not in error_dict['PIE']:
						error_dict['PIE'].append("PIE_TH")
				stepLst.append(int((adjustMax-adjust)/adjustStep+0.01))
				#cellsize steps
				#catch error, if cellrange+MinCellSize >= MaxCellSize
				if (cellRange >= MaxCellSize):
					cellRange = MinCellSize
					print "Max-MinCellSize was bigger than MaxCellSize -> set to MinCellSize to avoid unusable results! (" + str(MinCellSize) + ")
					error_dict['PIE'].append("PIE_CellSize")
					
				if cellRange < MinCellSize:
					cellRange = MinCellSize
					print "Max-MinCellSize was smaller than MinCellSize -> set to MinCellSize! (" + str(MinCellSize) + ")"
					if "PIE_CellSize" not in error_dict:
						error_dict['PIE'].append("PIE_CellSize")
	
				stepLst.append(int(floor((cellRange-MinCellSize)/cellStep)))
								
				stepLst = [x+1 for x in stepLst]
				#set step count to 1 if static values were chosen
				if minStatic:
					stepLst[0]=1
				if iteration_static:
					stepLst[1]=1
				if adjustStatic:
					stepLst[2]=1
				if conStatic:
					stepLst[3]=1
				
				print "Steps:", stepLst
				maxInterations = int(reduce(lambda x,y: x*y, stepLst))
				print "Iterations: ", maxInterations

				
				#start loop iterations	
				#loop order: 1. Circularity (i) 2. Iterations (j) 3. Threshold adjustment (k) 4. Contrast (l)
				#each loop alternates one of the not static values
				iteration = 1
				#use original values to get correct new values for each iteration
				#start_value_array.extend([MinCirc,iteration_steps,adjust,MinCellSize,CellDiameter])
				start_value_array[0], start_value_array[1], start_value_array[2], start_value_array[3] = MinCirc,iteration_steps,adjust,MinCellSize
				
				#Circularity
				for i in range(0,stepLst[0]):
					MinCirc = start_value_array[0] + MinCircStep * i 
					#Iterations
					for j in range(0,stepLst[1]):
						iteration_steps = start_value_array[1] + iteration_inc * j
						#Threshold adjustment
						for k in range(0,stepLst[2]):
							#get correct adjust value
							#calculate adjust value (adjust formula: if positive adjust = float(1 + adjust*0.01), if negative: adjust = float(float(1)/adjust)
							#this adjust holds the correct adjust for the current iteration
							this_adjust = start_value_array[2] + adjustStep * k
							adjust = float(float(1)+this_adjust*0.01)
							#Contrast
							for l in range(0,stepLst[3]):
								print "Iteration : ", iteration
								MinCellSize = start_value_array[3] + cellStep * l
								#save settings in settingsDictionary
								settingDic[iteration] = [MinCirc,iteration_steps,adjust,MinCellSize]
								#iteration algorithm
								for filename in os.listdir(folder):
									if filename.endswith(".JPG") or filename.endswith(".jpg") or filename.endswith(".TIFF") or filename.endswith(".tiff") or filename.endswith(".TIF") or filename.endswith(".tif"):
										
										imageValid = True
										currentImageName = filename
										#validate current image for statistical approach (xml file there?)
										imageValid = xml_reader(XMLfolder,filename)
										if imageValid:
											globalImageCount = globalImageCount + 1
											imp = IJ.openImage(os.path.join(folder,filename))
											#imp = WindowManager.getCurrentImage()
											# capture exception
											if imp is None:
												print "Could not open image from file:", filename
												if "FIE_loadImage" not in error_dict['FIE']:
													error_dict['FIE'].append("FIE_loadImage")
												continue		
											# get statistic values over function
											mean = getStatistics(imp)
											
										
											#dictionary saves mean values for every Image
											meanValues[filename] = round(mean,1)
											#CountDic saves for every Images the count
											CountDic[filename] = []
											
											###################################################
											pre_processed(imp)
											###################################################

								IJ.showProgress(iteration,maxInterations)
								#increase values after iteration steps
								iteration = iteration + 1
				globalImageCount = int(globalImageCount/maxInterations)
		#___________________________________ all other choices ___________________________________________			
			else:
				
				########## read single files and process them ############
				for filename in os.listdir(folder):
					if filename.endswith(".JPG") or filename.endswith(".jpg") or filename.endswith(".TIFF") or filename.endswith(".tiff") or filename.endswith(".TIF") or filename.endswith(".tif"):
						globalImageCount = globalImageCount + 1
						imageValid = True
						currentImageName = filename
						# Test if image is allowed to be processed (if no statistical values are required it's always processed)
					
						if imageValid:
							print "Processing", filename						
							imp = IJ.openImage(os.path.join(folder,filename))
							# capture exception
							if imp is None:
								print "Could not open image from file:", filename
								## log entry
								
								if "OE_loadImage" not in error_dict['OE']:
									error_dict['OE'].append("OE_loadImage")		
								continue		
							# get statistic values over function
							mean = getStatistics(imp)
							#dictionary saves mean values for every Image
							meanValues[filename] = round(mean,1)
							#CountDic saves for every Images the count
							CountDic[filename] = []
							
							###################################################
							pre_processed(imp)
							####################################################

	#___________________________________ Output functions ___________________________________________					
			if not(param_search) and len(CountDic):	
				
				b = datetime.now()
				finalTime = b-a	
				saveCountResults(CountDic)
			
			elif param_search and len(iterationDic) > 0:			
				
				if outlierChoice:
					#create brokenImageIteration Dictionary : saves for each iteration the broken images (names)
					brokenIterNameDic = {}
					
					brokenImageFilter()
					#use qualityScore again to get "corrected" values
					f1Score()								
				else:
					f1Score()
				
				b = datetime.now()
				finalTime = b-a
				#get input adjust value back
				#output
				saveTrainingResults()
			else:
				print "No images analysed! No output-file created!"
			print "Calculation finished... "
			#end calc time - for benchmarking
			
			#catch cancled dialogs
			if 'finalTime' in globals():
				print  str(finalTime.seconds) + "." + str(finalTime.microseconds/1000) + " Seconds needed to calculate all images."		
		#no dialog input	
		else: 
			print "No input from dialog"
	#if settingChoice is not None: 
	else:
		print "No input... Exiting program"
