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

baseDir='.'
inputMatrixFileName='contact.absCorr'
inputMatrixFilePath='/'.join([baseDir,inputMatrixFileName])

outputDir=baseDir
outputMatrixFileName='contact.CORR_BTW'
outputMatrixFilePath='/'.join([baseDir,outputMatrixFileName])

print "loading matrix data from %s"%(inputMatrixFilePath)
inputMatrixData=correlation_data_utilities.read_carma_matrix(
                        inputMatrixFilePath,has_header=True,returnHeader=True,verbose=False)
inputMatrix=correlation_data_utilities.corrDataDictToMat(inputMatrixData)

sources=[49]
targets=[432]
nResidues=inputMatrix.shape[0] #this should be a square matrix... if not we have problems


flowBtwMat=np.abs(correlation_data_utilities.getBtwMat(np.abs(inputMatrix),
                                                sources=sources,targets=targets,verbose=True))

correlation_data_utilities.write_mat_to_carma_matrix(
                        outputMatrixFilePath,flowBtwMat,verbose=True)

print 'done'
