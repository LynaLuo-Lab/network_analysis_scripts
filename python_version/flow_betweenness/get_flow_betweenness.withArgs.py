#run this file with python get_flow_betweenness.py

import numpy as np
import gc
import os
import sys
import subprocess
import itertools
import copy
import collections
import re
import glob
import pandas as pd

import correlation_data_utilities

if not (len(sys.argv)==5):
	print "usage: python get_flow_betweenness.withArgs.py inputMatrixFileName outputMatrixFileName sourceNodesString targetNodesString"
	print "sourceNodesString and targetNodesString should be quotted space separated lists of array indices, e.g. '1 2 3 4'"
	sys.exit()
else:
	baseDir='.'
	inputMatrixFileName=sys.argv[1]
	inputMatrixFilePath='/'.join([baseDir,inputMatrixFileName])
	
	outputDir=baseDir
	outputMatrixFileName=sys.argv[2]
	outputMatrixFilePath='/'.join([baseDir,outputMatrixFileName])
	
	print "loading matrix data from %s"%(inputMatrixFilePath)
	inputMatrixData=correlation_data_utilities.read_carma_matrix(
	                        inputMatrixFilePath,has_header=True,returnHeader=True,verbose=False)
	inputMatrix=correlation_data_utilities.corrDataDictToMat(inputMatrixData)
	
	sourcesString=sys.argv[3]
	sources=np.array(sourcesString.split(),dtype=int)
	targetsString=sys.argv[4]
	targets=np.array(targetsString.split(),dtype=int)
	nResidues=inputMatrix.shape[0] #this should be a square matrix... if not we have problems
	
	
	flowBtwMat=np.abs(correlation_data_utilities.getBtwMat(np.abs(inputMatrix),
	                                                sources=sources,targets=targets,verbose=True))
	
	correlation_data_utilities.write_mat_to_carma_matrix(
	                        outputMatrixFilePath,flowBtwMat,verbose=True)
	
	print 'done'
