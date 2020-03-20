from __future__ import print_function
import numpy as np
import pandas as pd
import scipy as sp
import re
import sys
import os
import collections
import copy
import tqdm
from tqdm import tqdm
from tqdm import tqdm_notebook
import networkx as nx 
from itertools import islice
from itertools import combinations
import gc

import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
import bokeh as bk
from bokeh.plotting import figure, output_notebook, output_file, reset_output, show, ColumnDataSource
from bokeh.models import LinearColorMapper, ColorBar, FuncTickFormatter, Title, BasicTicker
from bokeh.transform import linear_cmap, factor_cmap, LinearColorMapper
from bokeh.models.tools import HoverTool
from IPython.display import display, HTML

import nglview as nv

reset_output()
output_notebook()

noPytraj=False
try:
    import pytraj as pt
    from pytraj.utils.progress import ProgressBarTrajectory
except ImportError:
    print("could not import pytraj")
    noPytraj=True
    
import subprocess
import glob


def is_number(val):
    try:
        float(val)
    except ValueError:
        return False
    return True

def read_gCorrelation_data(filePath,verbose=False):
    with open(filePath) as inputFile:
        line = inputFile.readline()
        #lineTokens=line.split()
        while line:
            if("[" in line):
                foundStartLine=True
            if foundStartLine:
                lineTokens=line.split()
                if(lineTokens[0].isdigit() & lineTokens[2].isdigit()):
                    nRows=int(lineTokens[0])
                    nCols=int(lineTokens[2])
                break
        line=inputFile.readline()
        nEntries=nRows*nCols
        entries=np.zeros(nEntries)
        count=0
        if(nEntries > 1000):
            cMod=np.ceil(nEntries/100)
            lCount=0
        elif nEntries > 100:
            cMod=np.ceil(nEntries/10)
        else:
            cMod=1
        if verbose:
            print("Reading file "+filePath+":", end=' ')
        while line:
            lineTokens=line.split()
            lineEntries=np.extract(list(map(is_number,np.array(lineTokens))),np.array(lineTokens))
            entries[count:(count+len(lineEntries))]=lineEntries
            count=count+len(lineEntries)
            line=inputFile.readline()
            if verbose:
                if count % cMod == 0:
                    print(str(np.floor(1000 * count / nEntries)/10)+"%", end=' ')
        if verbose:
            print("")
    return collections.OrderedDict({"nRows":nRows,"nCols":nCols,"entries":entries})

def write_carma_matrix(filepath,nRows,nCols,dataEntries,
                       writeHeader=True,
                       header="((protein) and (not hydrogen)) and ((name CA) )",
                       verbose=False):
    if verbose:
        print("Writting carma matrix format file: "+filepath)
    with open(filepath,'w') as outputFile:
        if writeHeader:
            if verbose:
                print('Writting header')
            outputFile.write(header)
            outputFile.write('\n')
        if verbose:
            print('Writting '+str(nRows)+' rows:', end=' ')
            if nRows > 100:
                cMod=nRows/10
            elif nRows > 10000:
                cMod=nRows/100
                count=0
            else:
                cMod=1
        for iRow in np.arange(nRows):
            iStart=iRow*nCols
            iEnd=(iRow+1)*nCols
            outputFile.write(' '.join(map(str,dataEntries[iStart:iEnd])))
            outputFile.write('\n')
            if verbose and (iRow % cMod == 0):
                if cMod>1:
                    if (nRows > 10000) and (count > 9):
                        print("\,    ")
                        count=0
                    print("%4.1f%s "%((np.ceil((1000*iRow)/nRows)/10.0),"%"), end=' ')
                    if nRows > 10000:
                        count = count+1
                else:
                    print(".", end=' ')
        if verbose:
            print("\ndone")
            
            

def read_carma_matrix(filepath,has_header=True,returnHeader=False,verbose=False):
    #assumes all carma matrices are square matrices!
    if verbose:
        print("Reading carma matrix format file: "+filepath)
    with open(filepath,'r') as inputFile:
        if(has_header):
            if verbose:
                print('reading header')
            headerLine=inputFile.readline()
        line=inputFile.readline()
        count=0
        cMod=1
        while line:
            lineTokens=line.split()
            if(count==0):
                nRows=len(lineTokens)
                nEntries=nRows*nRows
                entries=np.zeros(nEntries)
                if(verbose):
                    print('Reading '+str(nRows)+' rows: ', end=' ')
                    if nRows > 1000:
                        cMod = np.ceil(nRows / 100)
                    elif nRows > 100:
                        cMod = np.ceil(nRows / 10)
                    else:
                        cMod = 1
            valArray=np.array(lineTokens)
            if (len(valArray)==nRows):
                iStart=count*nRows
                iEnd=iStart+len(valArray)
                sys.stdout.flush()
                entries[iStart:iEnd]=valArray
                if verbose & (count % cMod == 0):
                    if nRows > 10:
                        print(str(np.floor(1000*count/nRows)/10)+"% ", end=' ')
                    else:
                        print(".", end=' ')
                count=count+1
            line=inputFile.readline()
        if verbose:
            print("")
    if count > 0:
        if has_header & returnHeader:
            return collections.OrderedDict({"headerLine":headerLine,"nRows":nRows,"nCols":nRows,
                                            "entries":entries})
        else:
            return collections.OrderedDict({"nRows":nRows,"nCols":nRows,"entries":entries})
    else:
        print("Error! Data file appears empty!")

def convert_gCorrelationData_to_carmaMatrix(gCorrFilePath,outputFilePath,
                                            writeHeader=True,
                                            header="((protein) and (not hydrogen)) and ((name CA) )",
                                            verbose=False):
    #if verbose:
    #    print "Reading g_correlation data file"
    #    sys.stdout.flush()
    gCorrData = read_gCorrelation_data(filePath=gCorrFilePath,verbose=verbose)
    #if verbose:
    #    print "Writting g_correlation data to carma matrix format file"
    #    sys.stdout.flush()
    write_carma_matrix(filepath=outputFilePath,
                       nRows=gCorrData['nRows'],nCols=gCorrData['nCols'],
                       dataEntries=gCorrData['entries'],
                       writeHeader=writeHeader,header=header,verbose=verbose)
    if verbose:
        print("Conversion complete")

#Utilities for computing 'flow betweenness' scores for matrix
#representations of networks
def corrDataDictToMat(dataDict):
    return np.matrix(
        np.array(
            dataDict['entries']).reshape(dataDict['nRows'],dataDict['nCols']))

def corrMatToDataDict(mat):
    nRows=mat.shape[0]
    nCols=mat.shape[1]
    entries=np.array(mat).reshape(nRows*nCols)
    return(collections.OrderedDict({
                'nRows':nRows,
                'nCols':nCols,
                'entries':entries}))

def edgeDataDictToMatrix(dataDict,nRows=-1,nCols=-1):
    if nRows < 0:
        nRows=np.max(dataDict['Ei'])+1
    if nCols < 0:
        nCols=np.max(dataDict['Ej'])+1
    outMat=np.matrix(np.zeros([nRows,nCols])*0.0)
    outMat[dataDict['Ei'],dataDict['Ej']]=dataDict['Ew']
    return outMat

def write_dataDict_to_carma_matrix(filepath,dataDict,writeHeader=True,
                                       useDictHeader=True,
                                       header="((protein) and (not hydrogen)) and ((name CA) )",
                                       verbose=False):
    if not ( ('nRows' in dataDict) and \
             ('nCols' in dataDict) and \
             ('entries' in dataDict) ):
        print("ERROR! data dictionary is missing needed key value pairs. !!Aborting!!")
        return
    if (writeHeader and useDictHeader):
        if not ('headerLine' in dataDict) :
            print("WARNING! headerLine was missing from data dictionary.")
            print(" -Defaulting to: "+header)
            headerLine=header
        else:
            headerLine=dataDict['headerLine']
    write_carma_matrix(filepath=filepath,
                       nRows=dataDict['nRows'],nCols=dataDict['nCols'],
                       dataEntries=dataDict['entries'],writeHeader=writeHeader,
                       header=headerLine,verbose=verbose)
        
            

def write_mat_to_carma_matrix(filepath,mat,writeHeader=True,
                              header="((protein) and (not hydrogen)) and ((name CA) )",
                              verbose=False):
    dataDict=corrMatToDataDict(mat)
    if writeHeader:
        dataDict['headerLine']=header
    write_dataDict_to_carma_matrix(filepath,dataDict,
                                   writeHeader=writeHeader,useDictHeader=writeHeader,
                                   header=header,verbose=verbose)

#functions for generating network topology and edge weighting data with pytraj
    
def gen_perResCOM_traj(inputTrajPath,inputTopPath,residList,
                              resSelectMask='',resAtomMask='@CA',
                              COM_mask='',computeCommand='vector center :',
                              threads=1,verbose=False):
    #if noPytrajFlag:
    #    print "ERROR! gen_perResCOM_traj cannot be run if pytraj fails to load. Aborting"
    #    return
    #Default is to return per-residue center of mass mapped to alpha carbons of each residue
    if verbose:
        print("loading input trajectory")
    iterTraj=pt.iterload(inputTrajPath,top=inputTopPath)
    if verbose:
        print(iterTraj)
    commandList=["vector center :"+str(rid)+COM_mask for rid in residList]
    if verbose:
        if len(commandList) >= 10:
            print("first 10 mapped commands: ", end=' ')
            print(commandList[0:10])
        else:
            print("mapped commands: ", end=' ')
            print(commandList)
        print("running mapped commands")
    perResData=np.array(list(pt.compute(commandList,iterTraj,n_cores=threads).values()))
    if resSelectMask=='':
        residString=[str(rid) for rid in residList]
        resSelectMask=':'+','.join(residString)+resAtomMask
    if verbose:
        print("extracting residue selection subtrajectory")
    tempTraj=iterTraj[0:iterTraj.n_frames,resSelectMask]
    if verbose:
        print("subtrajectory info:")
        print(tempTraj)
    if not (tempTraj.shape[0]==perResData.shape[1] and \
            tempTraj.shape[1]==perResData.shape[0]):
        print("ERROR! mismatch between trajectory subselection coordinates and perResidue command data")
        print(" trajectory.xyz.shape=", end=' ')
        print(trajectory.xyz.shape, end=' ')
        print("; perResData.shape=", end=' ')
        print(perResData.shape[[1,0,2]])
        return
    if verbose:
        print("updating trajectory coordinates")
    for iDim in np.arange(3):
        tempTraj.xyz[:,:,iDim]=perResData[:,:,iDim].T
    print("done")
    return tempTraj

def gen_smoothed_contact_map(traj,resids,
                             timeAggFun=lambda timeSeries: 1.0*(np.mean(timeSeries)>.75),
                             distSmoothFun=lambda distArray:(6.0-np.clip(distArray,3.0,6.0))/(6.0-3.0),
                             verbose=False,verboseLevel=0):
    #if noPytrajFlag:
    #    print "ERROR! gen_perResCOM_traj cannot be run if pytraj fails to load. Aborting"
    #    return
    #traj should by a pytraj.trajectory or trajectory iterator
    #resids should be a list of residue ids (using AMBER numbering)
    #distSmoothFun should apply a smoothing kernel over a
    #len(resids) by traj.n_frames 2D np.ndarray
    #It defaults to a linear dropoff from an upper cutoff of 3.0 (yields 1.0 or full contact)
    #to a lower cutoff of 6.0 (yields 0.0 or no contact)
    #timeAggFun should apply an appropriate aggregation function
    #it defaults to computing the mean value and using a .75 high pass filter
    #to recover the unsmoothed contact map use distSmoothFun=lambda distArray: 1.0*(distArray<3.0)
    #verbose will turn on or off a simple text based progress monitor
    nResidues=len(resids)
    tempMat=np.zeros(nResidues)
    residArray=np.array(resids)
    if verbose:
        print("Computing contact matrix:", end=' ')
        if nResidues > 100:
            cMod=nResidues/100
        elif nResidues > 10:
            cMod=nResidues/10
        else:
            cMod=1
        lCount=0
    (mgrid_i,mgrid_j)=np.mgrid[0:nResidues,0:nResidues]
    tempMat=np.zeros([nResidues,nResidues])
    for iRow in np.arange(nResidues):
        if verbose:
            if((iRow % cMod == 0) or (iRow==(nResidues-1))):
                print("%4.1f%s "%(np.ceil(1000*iRow/nResidues)/10,"%"), end=' ')
                if lCount==10:
                    print("\n                         ", end=' ')
                    lCount=0
                else:
                    lCount=lCount+1
        rgrid_i=residArray[mgrid_i[iRow,:].flatten()]
        rgrid_j=residArray[mgrid_j[iRow,:].flatten()]
        commandList=list(map(lambda rid,rjd: 'nativecontacts :'+\
                        str(rid)+' :'+str(rjd)+' mindist',rgrid_i,rgrid_j))
        tempDists=np.array(list(pt.compute(commandList,traj).values())[2:(3*len(rgrid_i)):3])
        if verbose and (verboseLevel>0):
            print(commandList[0:3])
            print(type(tempDists), end=' ')
            print(" ", end=' ')
            print(tempDists.shape, end=' ')
            print(tempDists)
        matVals=[timeAggFun(timeSeries) for timeSeries in distSmoothFun(tempDists)]
        tempMat[mgrid_i[iRow,:].flatten(),mgrid_j[iRow,:].flatten()]=matVals
        tempMat[iRow,iRow]=0
    print("")
    return tempMat

#functions for generating current flow betweenness data from network matrices

def matLap(mat):
    if mat.shape[0]>mat.shape[1]:
        Lmat=copy.deepcopy(-mat.T).astype(float)
    else:
        Lmat=copy.deepcopy(-mat).astype(float)
    for iRow in np.arange(Lmat.shape[0]):
        Lmat[iRow,iRow]=0
        Lmat[iRow,iRow]=Lmat[iRow,iRow]-\
            np.sum(Lmat[iRow,:])
    if mat.shape[0]>mat.shape[1]:
        Lmat=Lmat.T
    return Lmat

def matAdj(mat):
    if mat.shape[0]>mat.shape[1]:
        Amat=copy.deepcopy(mat.T)
    else:
        Amat=copy.deepcopy(mat)
    Amat[np.arange(Amat.shape[0]),np.arange(Amat.shape[0])]=np.zeros(Amat.shape[0])
    #for iRow in np.arange(Amat.shape[0]):
    #    Amat[iRow,iRow]=0
    if mat.shape[0]>mat.shape[1]:
        Amat=Amat.T
    return Amat
    

def e_btw_from_Linv(Linv,Amat,sources,targets,verbose=False,verboseLevel=0,
                    useProgressBar=True,pbarFun=tqdm):
    #This algorithm makes use of the full moore-penrose pseudo-inverse
    #which is generally quite dense. The algorithm will scale quite poorly
    #for very large systems.
    #The bottleneck step (after computation of the pseudo inverse) scales as O(m*s*t) where: 
    #   m is the number of non-zero ntries in the adjacency matrix (network edges)
    #   s is the number of source nodes
    #   t is the number of target nodes
    #
    #This could yield O(n^4) for networks of n nodes in the worst case!
    #
    #This is only correct for the case where sources and targets are disjoint
    #sets. If not, the scaling factor must be adjusted to be equal the number of
    #unique combinations of sources and targets.
    eMat=copy.deepcopy(Amat).astype(float)
    #some basic sanity checks
    if((Linv.shape[0]!=Linv.shape[1]) or 
       (Amat.shape[0]!=Amat.shape[1])):
        print("ERROR! Input matrices must by square!")
        eMat[:,:]=0
        return(eMat)
    if((Linv.shape[0]!=Amat.shape[0]) or
       (Linv.shape[1]!=Amat.shape[1])):
        print("ERROR! Input matrices must have the same shape!")
        eMat[:,:]=0
        return(eMat)
    if ((np.min(sources)<0) or
        (np.min(targets)<0) or
        (np.max(sources)>=Linv.shape[0]) or
        (np.max(targets)>=Linv.shape[1])):
        print("ERROR! invalid source or target index detected!")
        eMat[:,:]=0
        return(eMat)
    #get indices of edges... e.g. non-zero entries of Amat
    (Ei,Ej)=np.nonzero(Amat)
    if verbose:
        if verboseLevel > 2:
            print("Linv:", end=' ')
            print(Linv)
            print("Amat:", end=' ')
            print(Amat)
        print("computing betweenness for %g edges"%len(Ei))
        if verboseLevel > 0:
            print("Ei:", end=' ')
            print(Ei)
            print("Ej:", end=' ')
            print(Ej)
    if verbose and useProgressBar:
        Ebtw=np.array(list(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum([np.sum([np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]) for trg in targets]) for src in sources]),
                          pbarFun(Ei,file=sys.stdout),
                          Ej)))/(len(sources)*len(targets))
    else:
        Ebtw=np.array(list(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum([np.sum([np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]) for trg in targets]) for src in sources]),
                          Ei,
                          Ej)))/(len(sources)*len(targets))
    if verbose:
        if verboseLevel > 0:
            print("Ebtw:", end=' ')
            print(Ebtw)
            if verboseLevel > 1:
                print("(Ei,Ej):Ebtw;eMat")
    for iInd in np.arange(len(Ebtw)):
        eMat[Ei[iInd],Ej[iInd]]=Ebtw[iInd]
        if verbose :
            if verboseLevel > 1:
                print("(", end=' ')
                print(Ei[iInd], end=' ')
                print(",", end=' ')
                print(Ej[iInd], end=' ')
                print("):", end=' ')
                print(Ebtw[iInd], end=' ')
                print(";", end=' ')
                print(eMat[Ei[iInd],Ej[iInd]])
    return(eMat)
    

def getBtwMat(mat,sources,targets,verbose=False,verboseLevel=0,
              useProgressBar=False,pbarFun=tqdm):
    #Given a (possibly weighted) network in matrix format (mat)
    #and a set of source and target nodes (sources and targets)
    #return the corresponding network with flow betweenness
    #edge weights.
    #The pseudo inverse step is the most likely bottleneck, however
    #the betweenness calculation iteself scales as O(m*s*t)
    #where m=number of network edges, s=number of source nodes, and
    #t=number of target nodes. 
    #At worst case, this could yield O(n^4) for an n-node matrix!
    #Also, the sources and targets must be disjoint sets or the results
    #will be incorrect.
    if verbose:
        print("computing matrix Laplacian")
    Lmat=matLap(copy.deepcopy(mat))
    if verbose:
        print("extracting weighted adjacency matrix")
    Amat=matAdj(copy.deepcopy(mat))
    if verbose:
        print("computing moore-penrose inverse of matrix Laplacian")
    Linv=np.linalg.pinv(Lmat)
    if verbose:
        print("generating flow betweenness scores")
    return(e_btw_from_Linv(Linv,Amat,sources,targets,
                           verbose=verbose,verboseLevel=verboseLevel,
                           useProgressBar=useProgressBar,
                           pbarFun=pbarFun))
    
    

def calcCorrDissipation(corrMat,btwMat):
    return np.sum(np.abs(np.array(btwMat))*\
                  np.abs(np.array(btwMat))/np.abs(np.array(corrMat)))

#functions for running and parsing subopt files generated via the VMD 'subopt' command

def validate_subopt_file(filepath,verbose=False):
    if not os.path.exists(filepath):
        if verbose:
            print(filepath+" does not exist.")
        return False
    if not os.path.isfile(filepath):
        if verbose:
            print(filepath+" is not a file.")
        return False
    foundPathStart=False
    foundPathCount=False
    with open(filepath,'r') as suboptFile:
        for line in suboptFile:
            if re.search('The final paths are',line):
                foundPathStart=True
            if re.search('Number of paths is',line):
                foundPathCount=True
            if foundPathStart and foundPathCount:
                break
    if foundPathStart and foundPathCount:
        if verbose:
            print(filepath+' is a valid subopt file')
        return True
    else:
        if verbose:
            print(filepath+" : ", end=' ')
            if not foundPathStart:
                print("is missing paths section, ", end=' ')
            if not foundPathCount:
                print("is missing path count")
        return False

def get_subopt_pathCount(filepath,verbose=False):
    if not validate_subopt_file(filepath,verbose=verbose):
        print("ERROR! "+filepath+" is not a valid subopt file.")
        return -1
    else:
        with open(filepath,'r') as suboptFile:
            for line in suboptFile:
                if re.search('Number of paths is',line):
                    tokens=str.split(line)
                    pathCount=tokens[len(tokens)-1]
                    if verbose:
                        print('path count = '+str(pathCount))
                    if not str.isdigit(pathCount):
                        break
                    return pathCount
        print("ERROR! Something went wrong, the file seemed valid but had an invalid path count line")
        return -1

def get_subopt_pathData(filepath,verbose=False):
    pathData=collections.OrderedDict()
    if not validate_subopt_file(filepath,verbose=verbose):
        print("ERROR! "+filepath+" is not a valid subopt file.")
        return pathData
    else:
        foundPathStart=False
        foundPathCount=False
        pathData['paths']=[]
        pathData['lengths']=[]
        with open(filepath,'r') as suboptFile:
            for line in suboptFile:
                if re.search('Number of paths is',line):
                    foundPathCount=True
                    tokens=str.split(line)
                    pathCount=tokens[len(tokens)-1]
                    if verbose:
                        print('path count = '+str(pathCount))
                    if not str.isdigit(pathCount):
                        foundPathCount=False
                        pathCount=-1
                    break
                if not foundPathStart:
                    if re.search('The final paths are',line):
                        foundPathStart=True
                else:
                    tokens=list(map(int,re.sub('[,)(]','',line).split()))
                    tempPath=tokens[0:(len(tokens)-1)]
                    pathData['paths'].append(tempPath)
                    pathData['lengths'].append(tokens[len(tokens)-1])
        if not foundPathStart:
            print("ERROR! "+filepath+" seemed valid but path section was apparently absent")
            return pathData
        if not foundPathCount:
            print("Warning! final path count line was missing or ill-formed!")
            pathData['count']=len(pathData['paths'])
        else:
            if len(pathData['paths']) != int(pathCount):
                print("Warning! subopt file lists number of paths as", end=' ')
                print(str(pathCount), end=' ')
                print("but", end=' ')
                print(len(pathData['paths']), end=' ')
                print("paths were found.")
                pathData['count']=len(pathData['paths'])
            else:
                pathData['count']=pathCount
        return pathData

def run_external_subopt(networkMatFilePath,outputFilePath,
                        sourceNode,targetNode,dilationValue,
                        externalSuboptCommand='subopt',
                        returnSuboptData=False,
                        verbose=False):
    if sourceNode > targetNode:
        sID=targetNode
        tID=sourceNode
    else:
        sID=sourceNode
        tID=targetNode
    suboptCommand=' '.join([
            externalSuboptCommand,
            networkMatFilePath,
            outputFilePath,
            str(dilationValue),str(sID),str(tID)])
    if verbose:
        print('running external subopt command: '+suboptCommand)
    os.system(suboptCommand)
    if returnSuboptData:
        return get_subopt_pathData(outputFilePath+".out",verbose=verbose)
    

def run_subopt_till_pathCount(networkMatFilePath,sourceNode,targetNode,
                              outputDir,outputBaseName='subopt',
                              minPathCount=0,percentDilationIncrement=10,
                              externalSuboptCommand='subopt',
                              returnSuboptData=False,onlyFinalRun=True,
                              verbose=False,verboseLevel=0):
    subVerbose=(verbose and (verboseLevel > 0))
    if verbose:
        print('running iterative dilation till '+str(minPathCount)+' paths are attained:')
        if not subVerbose:
            print('(%dilation,pathCount):', end=' ')
    percentDilation=0
    dilationValue=0
    outputFileName='.'.join([
            outputBaseName,
            '_'.join([str(sourceNode),str(targetNode)]),
            '_'.join(['dilation',str(percentDilation)])])
    outputFilePath='/'.join([outputDir,outputFileName])
    run_external_subopt(networkMatFilePath,outputFilePath,
                        sourceNode,targetNode,dilationValue,
                        externalSuboptCommand,returnSuboptData=False,
                        verbose=subVerbose)
    suboptDataFilePath='.'.join([outputFilePath,'out'])
    if not validate_subopt_file(suboptDataFilePath):
        print("ERROR! subopt failed to generate a valid output file. Aborting")
        if returnSuboptData:
            if onlyFinalRun:
                return collections.OrderedDict()
            else :
                return [collections.OrderedDict()]
        else:
            return
    pathCount=get_subopt_pathCount(suboptDataFilePath)
    if verbose:
        print("(%g,%g)"%(float(percentDilation),float(pathCount)), end=' ')
        if not subVerbose:
            print(",", end=' ')
        else:
            print("")
    tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)
    optPathLength=np.min(list(map(float,tempData['lengths'])))
    if returnSuboptData:
        if onlyFinalRun:
            suboptData=tempData
        else:
            suboptData=[tempData]
    while float(pathCount) < float(minPathCount):
        percentDilation=percentDilation+percentDilationIncrement
        dilationFactor=percentDilation/100.0
        dilationValue=optPathLength*dilationFactor
        outputFileName='.'.join([
            outputBaseName,
            '_'.join([str(sourceNode),str(targetNode)]),
            '_'.join(['dilation',str(percentDilation)])])
        outputFilePath='/'.join([outputDir,outputFileName])
        run_external_subopt(networkMatFilePath,outputFilePath,
                            sourceNode,targetNode,dilationValue,
                            externalSuboptCommand,returnSuboptData=False,
                            verbose=subVerbose)
        suboptDataFilePath='.'.join([outputFilePath,'out'])
        if not validate_subopt_file(suboptDataFilePath):
            print("ERROR! subopt failed to generate a valid output file. Aborting")
            if returnSuboptData:
                return suboptoptData
            else:
                return
        pathCount=get_subopt_pathCount(suboptDataFilePath)
        if verbose:
            print("(%g,%g)"%(float(percentDilation),float(pathCount)), end=' ')
        if not subVerbose:
            print(",", end=' ')
        else:
            print("")
        tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)   
        if returnSuboptData:
            if onlyFinalRun:
                suboptData=tempData
            else:
                suboptData.append(tempData)
    if verbose and not subVerbose:
        print("")
    if verbose:
        print("DONE!")
    if returnSuboptData:
        return suboptData
    
def get_subopt_dilations_data(suboptDir,basePattern,sep='.',
                              onlyMaxDilation=True,includeDilationValues=False,
                              includeFileNames=False,
                              verbose=False,verboseLevel=0):
    fileSearchPattern=str(sep).join([basePattern,'_'.join(['dilation','*']),'out'])
    if verbose:
        print('file search pattern = '+fileSearchPattern)
    searchPathPattern='/'.join([suboptDir,fileSearchPattern])
    if verbose:
        print('searchPathPattern = '+searchPathPattern)
    filePathList=glob.glob(searchPathPattern)
    fileNameList=[filepath.split('/')[-1] for filepath in filePathList]
    if verbose and (verboseLevel > 0):
        print('file name list: '+'\n'.join(fileNameList))
    if not onlyMaxDilation:
        suboptDataSets=[]
    else:
        suboptData=collections.OrderedDict()
        maxDilationValue=-1
    for iFile in np.arange(len(filePathList)):
        suboptFilePath=filePathList[iFile]
        suboptFileName=fileNameList[iFile]
        if validate_subopt_file(suboptFilePath):
            suboptData=get_subopt_pathData(suboptFilePath)
            nameTokens=suboptFileName.split(sep)
            dilationToken=[token for token in nameTokens if 'dilation' in token][-1]
            dilationValue=int(dilationToken.split("_")[-1])
            if includeDilationValues:
                suboptData['dilation']=dilationValue
            if includeFileNames:
                suboptData['filename']=suboptFileName
            if onlyMaxDilation:
                if dilationValue > maxDilationValue:
                    suboptDataSets=suboptData
                    maxDilationValue=dilationValue
            else:
                suboptDataSets.append(suboptData)
        else:
            print("Warning: "+suboptFileName+" was not a valid subopt file.")
    if onlyMaxDilation:
        return suboptData
    else:
        return suboptDataSets

def get_index_of_nth_maxRank_element(valList,n,verbose=False):
    #this assumes that valList is already sorted in ascending order!
    u,v=np.unique(valList,return_inverse=True)
    maxRanks=(np.cumsum(np.bincount(v,minlength=u.size))-1)[v]
    if verbose:
        print('value array:', end=' ')
        print(valList)
        print('rank array:', end=' ')
        print(maxRanks)
    if np.max(maxRanks) < n:
        print("get nth maxRank: Warning! there are not enough elements")
        return(len(valList)-1)
    else:
        return [m for m,i in enumerate(maxRanks) if i >= n][0]

def get_top_n_pathData_paths(pathData,n,verbose=False,verboseLevel=0):
    subVerbose=(verbose and (verboseLevel>0))
    outData=copy.deepcopy(pathData)
    pathSortingArray=np.argsort(outData['lengths'])
    outData['paths']=[outData['paths'][i] for i in pathSortingArray]
    outData['lengths']=[outData['lengths'][i] for i in pathSortingArray]
    maxPathIndex=get_index_of_nth_maxRank_element(outData['lengths'],n,
                                                  verbose=subVerbose)
    if verbose:
        print('max path index = '+str(maxPathIndex))
    outData['paths']=outData['paths'][0:(maxPathIndex+1)]
    outData['lengths']=outData['lengths'][0:(maxPathIndex+1)]
    outData['count']=len(outData['lengths'])
    return outData

def get_pathData_node_count_array(pathData,nNodes,inputIndexBase=0,verbose=False):
    outArray=np.zeros(nNodes)
    for path in pathData['paths']:
        outArray[path-inputIndexBase]=outArray[path-inputIndexBase]+1.0
    return outArray

def get_pathData_node_frequency_array(pathData,nNodes,inputIndexBase=0,verbose=False):
    outArray=get_pathData_node_count_array(pathData,nNodes,inputIndexBase,verbose)/\
        pathData['count']
    return outArray
    

def get_pathData_edge_count_matrix(pathData,nNodes,inputIndexBase=0,verbose=False):
    outMat=np.matrix(np.zeros([nNodes,nNodes]))
    for path in pathData['paths']:
        outMat[path[0:(len(path)-1)]-inputIndexBase,path[1:len(path)]-inputIndexBase]=outMat[
            path[0:(len(path)-1)]-inputIndexBase,path[1:len(path)]-inputIndexBase]+1.0
        outMat[path[1:(len(path))]-inputIndexBase,path[0:(len(path)-1)]-inputIndexBase]=outMat[
            path[1:(len(path))]-inputIndexBase,path[0:(len(path)-1)]-inputIndexBase]+1.0
    return outMat

def get_pathData_edge_frequency_matrix(pathData,nNodes,inputIndexBase=0,verbose=False):
    outMat=get_pathData_edge_count_matrix(pathData,nNodes,inputIndexBase,verbose)/\
        pathData['count']
    return outMat

def serialize_pathData_lengths(pathData):
    return ','.join(map(str,pathData['lengths']))

def serialize_pathData_paths(pathData):
    return ','.join(
        ['_'.join(map(str,pathArray)) for pathArray in pathData['paths']])

def get_1Darray_maxRanks(valArray,unsorted=True,invert=False,verbose=False):
    if unsorted:
        sorting_Array=np.argsort(valArray)
        sortedVals=[valArray[sInd] for sInd in sorting_Array]
        desorting_Array=np.zeros(len(sorting_Array),dtype=int)
        desorting_Array[sorting_Array]=np.arange(len(sorting_Array))
    else:
        sortedVals=valArray
    u,v=np.unique(sortedVals,return_inverse=True)
    if verbose:
        print('u:', end=' ')
        print(u)
        print('v:', end=' ')
        print(v)
    maxRanks=(np.cumsum(np.bincount(v,minlength=u.size))-1)[v]
    if unsorted:
        if invert:
            maxRank=np.max(maxRanks)
            outVals=[maxRank-maxRanks[sInd] for sInd in desorting_Array]
        else:
            outVals=[maxRanks[sInd] for sInd in desorting_Array]
        return outVals
    else:
        if invert:
            maxRank=np.max(maxRanks)
            return [maxRank-maxRanks[sInd] for sInd in np.arange(len(maxRanks))]
        else:
            return maxRanks
        
def get_matrix_element_maxRankings(mat,invert=False,verbose=False):
    unfoldedMat=np.array(copy.deepcopy(mat)).reshape(np.prod(mat.shape))
    unfoldedRankingMat=np.array(get_1Darray_maxRanks(
            unfoldedMat,invert=invert,verbose=verbose))
    rankingMat=unfoldedRankingMat.reshape(mat.shape)
    return rankingMat

#functions to facilitate loading path data from WISP logfiles
def load_wispLog_paths(logFilePath):
    pathList=[]
    with open(logFilePath,'r') as logFile:
        foundPathStart=False
        line=logFile.readline()
        while line and not foundPathStart:
            if 'Output identified paths' in line:
                foundPathStart=True
                #print line
            line=logFile.readline()
        if not foundPathStart:
            sys.stderr.write("Error: End of log file reached before path section was found")
        else:
            iPath=-1
            while '#' in line:
                if 'Path' in line:
                    if iPath>=0:
                        pathList.append(tempPath)
                        tempPath=[]
                        iPath=iPath+1
                    else:
                        tempPath=[]
                        iPath=iPath+1
                elif ('Length' in line):
                    tempPath.append(line.split()[2])
                elif ('Nodes:' in line):
                    tokens=np.array(line.split(),dtype='|S')
                    for node in tokens[2:len(tokens):2]:
                        tempPath.append(node)
                line=logFile.readline()
    return pathList

def simple_wispNode_converter(wispNode):
    return int(wispNode.split('_')[len(wispNode.split('_'))-1])

def wispPath_to_nodeIndPath(wispPath,
                            wispNode_to_nodeInd_function=simple_wispNode_converter):
    nodeIndPath=np.zeros(len(wispPath)-1)
    for iNode,wispNode in enumerate(wispPath[1:]):
        nodeIndPath[iNode]=wispNode_to_nodeInd_function(wispNode)
    return np.array(nodeIndPath,dtype=int)

def wispPaths_to_pathData(wispPaths,
                          wispNode_to_nodeInd_function=simple_wispNode_converter):
    pathData={
        'lengths':[],
        'paths':[],
        'count':len(wispPaths)
    }
    for iPath,wispPath in enumerate(wispPaths):
        pathData['lengths'].append(float(wispPath[0]))
        pathData['paths'].append(
            wispPath_to_nodeIndPath(wispPath,
                                    wispNode_to_nodeInd_function))
    return pathData

def load_wispLog_pathData(logFilePath,
                          wispNode_to_nodeInd_function=simple_wispNode_converter):
    return wispPaths_to_pathData(
        load_wispLog_paths(logFilePath),
        wispNode_to_nodeInd_function)

def netMatDict_to_nodeDataTable(netMatDict,
                                keySep='.',keyColNames=None,
                                indexCols=None,nodeIndName='X',
                                nodeInds=None,lapFunc=(lambda x: 2)):
    nodeTables=[]
    if (keyColNames is None):
        keyNames=['Key_%g'%iPart for iPart,part in \
                  enumerate(list(netMatDict.keys())[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
    with tqdm(len(netMatDict)) as pbar:
        for matKey in list(netMatDict.keys()):
            pbar.set_description_str(matKey)
            keyParts=matKey.split(keySep)
            tempMatLap=matLap(netMatDict[matKey])
            nodeVals=np.array(np.matrix(tempMatLap).diagonal()).flatten()
            norm=np.array(
                [lapFunc(
                    np.concatenate([tempMatLap[iRow,0:(iRow-1)],
                                    tempMatLap[iRow,(iRow+1):]])
                 ) for iRow,row in enumerate(tempMatLap)])
            nodeVals=nodeVals/norm
            nNodes=len(nodeVals)
            if nodeInds is None:
                inds=np.arange(nNodes)
            tempTable=pd.DataFrame({
                nodeIndName:inds,
                'value':nodeVals
            })
            for iCol,keyCol in enumerate(keyNames):
                tempTable[keyCol]=[keyParts[iCol]]*nNodes
            nodeTables.append(copy.deepcopy(tempTable))
            pbar.update()
    nodeTable=pd.concat(nodeTables)
    nodeDataTable=nodeTable.pivot_table(
        index=indCols,columns=nodeIndName,
        values='value',aggfunc=np.mean,fill_value=0
    )
    nodeDataTable.columns=np.array(nodeDataTable.columns)
    nodeDataTable=nodeDataTable.reset_index()
    return nodeDataTable

def netMatDict_to_edgeDataTable(netMatDict,
                                keySep='.',keyColNames=None,
                                indexCols=None,edgeIndNames=['X','Y'],
                                nodeInds=None,sparse=True):
    edgeTables=[]
    if (keyColNames is None):
        keyNames=['Key_%g'%iPart for iPart,part in \
                  enumerate(list(netMatDict.keys())[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
    with tqdm(len(netMatDict)) as pbar:
        for matKey in list(netMatDict.keys()):
            pbar.set_description_str('%s'%matKey)
            keyParts=matKey.split('.')
            tempMat=netMatDict[matKey]
            if sparse:
                edgeInds=np.nonzero(tempMat)
            else:
                pairs=np.array([[ii,jj] \
                       for ii in np.arange(tempMat.shape[0]) \
                       for jj in np.arange(tempMat.shape[1])])
                edgeInds=(pairs[:,0],pairs[:,1])
            edgeVals=np.array(tempMat)[edgeInds]
            nEdges=len(edgeVals)
            tempTable=pd.DataFrame({
                edgeIndNames[0]:edgeInds[0],
                edgeIndNames[1]:edgeInds[1],
                'value':edgeVals
            })
            for iCol,keyCol in enumerate(keyNames):
                tempTable[keyCol]=[keyParts[iCol]]*nEdges
            edgeTables.append(copy.deepcopy(tempTable))
            pbar.update()
    edgeTable=pd.concat(edgeTables)
    edgeDataTable=edgeTable.pivot_table(
        index=indCols,columns=edgeIndNames,
        values='value',aggfunc=np.mean,fill_value=0)
    edgeDataTable.columns=edgeDataTable.columns.map(
        lambda x: '_'.join(['%g'%xi for xi in x]))
    edgeDataTable=edgeDataTable.reset_index()
    return edgeDataTable

#taken directly from the networkx manual
def k_shortest_paths(G, source, target, k, weight=None):
     return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))
    
def converge_subopt_paths_betweenness(inputNetwork,source,target,weight='weight',
                                maxPaths=100,tolerance=1e-6,giveAlphas=False,verbose=False):
    '''Take additional paths between a source / target pair until the betweenness
       centrality of nodes within those paths computed over the attained suboptimal paths
       converges to a given tolerance.
       the convergence critera 'alpha' is computed as:
       alpha = (sum(current path node usage counts) / (number of paths + sum(all path node usage counts))
       the function will return a list of all generated paths. If the "giveAlphas" option is
       turned on, it will also return a list of the alpha after each iteration (path generated)
       this is useful for comparing convergence over different weighting schemes.
       if the "verbose" option is turned on, the alpha after each iteration will be printed
       to standard out as calculation proceeds.'''
    nNodes=len(inputNetwork.nodes())
    pathGenerator=nx.shortest_simple_paths(inputNetwork, source, target, weight=weight)
    nodeCounts=np.zeros(nNodes)
    pathList=[]
    alphas=[]
    newPath=next(pathGenerator)
    pathList.append(copy.deepcopy(newPath))
    pathCounts=np.unique(newPath,return_counts=True)
    tempCounts=np.zeros(nNodes)
    tempCounts[pathCounts[0]]=pathCounts[1]
    nPaths=1
    alpha=np.sum(tempCounts)/np.sum(nodeCounts+tempCounts)/nPaths*1.
    alphas.append(alpha)
    if verbose:
            print("%.3e"%alpha,end=", ")
    while (alpha>tolerance) & (nPaths<maxPaths):
        nodeCounts+=tempCounts
        
        newPath=next(pathGenerator)
        pathList.append(copy.deepcopy(newPath))
        pathCounts=np.unique(newPath,return_counts=True)
        tempCounts=np.zeros(nNodes)
        tempCounts[pathCounts[0]]=pathCounts[1]
        nPaths+=1
        alpha=np.sum(tempCounts)/np.sum(nodeCounts+tempCounts)/nPaths*1.0
        alphas.append(alpha)
        if verbose:
            print("%.3e"%alpha,end=", ")
    if alpha>tolerance:
        print("Maximum number of paths reached before converging betweenness score")
        print("Last relative count delta: %.3e"%alpha)
    if giveAlphas:
        return((pathList,alphas))
    else:
        return(pathList)
    
#A convenience function for calculating the length of an arbitrary path
#in a weighted graph
def calculatePathLength(pathGraph,path,weight='weight'):
    return(np.sum([pathGraph.edges()[(edge[0],edge[1])][weight] \
                   for edge in zip(path[:-1],path[1:])]))

#Utilities for computing distance topology in pytraj
def checkCollision(Rmin1,Rmax1,Rmin2,Rmax2,collisionRadius=0,axis=None):
    return(np.product((Rmin1-collisionRadius)<=(Rmax2+collisionRadius),axis=axis) * \
           np.product((Rmax1+collisionRadius)>=(Rmin2-collisionRadius),axis=axis))

def collisionCount(Rmin1,Rmax1,Rmin2,Rmax2,
                   collisionRadius=0,crdAxis=1):
    return(np.sum(np.product((Rmin1-collisionRadius)<=(Rmax2+collisionRadius),axis=crdAxis) * \
                  np.product((Rmax1+collisionRadius)>=(Rmin2-collisionRadius),axis=crdAxis)))

#This just returns whether or not a collision has ever happened, only slightly faster
#than the series apparently, though likely smaller memory requirements
def compute_BoxCollision_matrix(traj,collisionRadius=0.,resInds=None,showProgress=False):
    if resInds is None:
        nRes=traj.top.n_residues()
        resnums=np.arange(nRes)
    else:
        resnums=resInds
        nRes=len(resInds)
    
    resInds=[traj.topology.atom_indices(':%g'%iRes) for iRes in resnums+1]
    
    if showProgress:
        resIter=tqdm_notebook(resInds,desc='Computing Residue Minimum Bounds')
    else:
        resIter=resInds
    resMinBounds=np.array([np.min(traj.xyz[:,resInd,:],axis=1) \
                           for resInd in resIter])
    
    if showProgress:
        resIter=tqdm_notebook(resInds,desc='Computing Residue Maximum Bounds')
    else:
        resIter=resInds
    resMaxBounds=np.array([np.max(traj.xyz[:,resInd,:],axis=1) \
                           for resInd in resIter])
    
    resPairs=np.array(list(combinations(np.arange(nRes),2)))
    if showProgress:
        pairIter=tqdm_notebook(resPairs,desc='Computing Collisions')
    else:
        pairIter=resPairs
    collisionCheckArray=[
        checkCollision(resMinBounds[resPair[0]],resMaxBounds[resPair[0]],
                       resMinBounds[resPair[1]],resMaxBounds[resPair[1]],
                       collisionRadius=collisionRadius) \
        for resPair in pairIter]

    collisionMat=sp.sparse.coo_matrix(
        (collisionCheckArray,
         (resPairs[:,0],resPairs[:,1])),shape=(nRes,nRes))
    collisionMat=collisionMat+collisionMat.T
    return(collisionMat)

#Counts the number of frames where each residue pair has collided
#returns the result as a sparse matrix (scipy coo format)
def compute_BoxCollision_CountMatrix(traj,collisionRadius=0.,
                                     resinds=None,
                                     minBounds=None,maxBounds=None,
                                     showProgress=False,
                                     frameAxis=0,indAxis=1,crdAxis=2,
                                     returnBoundVecs=False):
    if resinds is None:
        nRes=traj.top.n_residues
        resnums=np.arange(nRes)
        resInds=[traj.topology.atom_indices(':%g'%iRes) for iRes in resnums+1]
    else:
        #resnums=resInds
        nRes=len(resinds)
        resInds=resinds
    
    #print(len(resInds))
    if (minBounds is None): #| (len(minBounds) != len(resInds)):
        if showProgress:
            resIter=tqdm_notebook(resInds,desc='Computing Residue Minimum Bounds')
        else:
            resIter=resInds
        resMinBounds=np.array([np.min(traj.xyz[:,resIndSet,:],axis=indAxis) \
                               for resIndSet in resIter])
    else:
        resMinBounds=minBounds
    
    if (maxBounds is None): #| (len(maxBounds) != len(resInds)):
        if showProgress:
            resIter=tqdm_notebook(resInds,desc='Computing Residue Maximum Bounds')
        else:
            resIter=resInds
        resMaxBounds=np.array([np.max(traj.xyz[:,resIndSet,:],axis=indAxis) \
                               for resIndSet in resIter])
    else:
        resMaxBounds=maxBounds
    
    resPairs=np.array(list(combinations(np.arange(nRes),2)))
    if showProgress:
        pairIter=tqdm_notebook(resPairs,desc='Computing Collisions')
    else:
        pairIter=resPairs
    #print(resMinBounds.shape)
    collisionCheckArray=[
        collisionCount(resMinBounds[resPair[0]],resMaxBounds[resPair[0]],
                       resMinBounds[resPair[1]],resMaxBounds[resPair[1]],
                       collisionRadius=collisionRadius,crdAxis=crdAxis-1) \
        for resPair in pairIter]

    #return(collisionCheckArray)
    collisionMat=sp.sparse.coo_matrix(
        (np.concatenate([collisionCheckArray,collisionCheckArray]),
         (np.concatenate([resPairs[:,0],resPairs[:,1]]),
          np.concatenate([resPairs[:,1],resPairs[:,0]]))),shape=(nRes,nRes))
    if returnBoundVecs:
        return(collisionMat,resMinBounds,resMaxBounds)
    else:
        return(collisionMat)
    
def compute_pairwise_minDist_data(traj,
                                  resindexpairs=None,
                                  chunkSize=1000,
                                  outFilePath=None,
                                  returnData=True,
                                  showProgress=False):
    '''
        traj: a pytraj trajectory
        resindexpairs: pairs of residues to compute distance series for. These should be base 0 indices.
            The default will compute all possible residue pair combinations. This option can allow you
            to filter which pairs to compute. E.g. if you have an nResidue by nResidue matrix (filterMatrix)
            wich has nonzero values for only the pairs to be computed you could use:
            resindexpairs=zip(np.nonzero(filterMatrix)[0],np.nonzero(filterMatrix)[1])
        showProgress: if set to True, display a tqdm_notebook style progress bar
        chunkSize: Due to technical considerations, it is apparently faster to do several hundred to
            several thousand pair distance calculations at a time. This controls the number of
            pair distances calculated in one pass over the entire trajectory. For the IGPS system
            1000 seemed to be roughly optimal, but the optimal value will likely vary depending on system size
            and trajectory length.
        outFilePath: if set, the data will be written chunk by chunk to the specified filepath as it is
            generated.
        returnData: return an array containing the computed data. If turned off, nothing will be returned.
            This can be useful when handling a very large number of pairs (i.e. if the computed data
            will not fit in memory). Setting this to false and providing an outFilePath will cause the 
            data to be written directly to disk instead of stored in an in-memory array.
            Note that you will still need to be able to fit at least pair set chunk worth of data in memory.
        
        returns: an M-Residue_Pair by N-Frame array where rows are the residue pair and columns are 
            trajectory frames. Each entry is the minimum residue-residue interatomic distance for the
            given residue pair at the given frame.
    '''
    if resindexpairs is None:
        resIndices=np.arange(traj.Top.n_residues)
        resIndexPairs=list(combinations(resIndices,2))
    else:
        resIndices=resindexpairs
    
    distData=[]
    
    chunkStart=0
    nPairs=len(resIndices)
    
    if showProgress:
        pbar=tqdm_notebook(total=nPairs,
                                desc='computing pair distances')
    count=0
    while chunkStart<nPairs:
        chunkEnd=chunkStart+chunkSize
        pairIter=resIndices[chunkStart:chunkEnd]
        commandList=['nativecontacts mindist :{:g} :{:g}'.format(resPair[0]+1,resPair[1]+1) \
                     for resPair in pairIter]
        tempData=list(
            pt.compute(commandList,traj).values())[2::3]
        if returnData:
            distData.append(copy.deepcopy(tempData))
        chunkStart=chunkEnd
        if not (outFilePath is None):
            pbar.set_description('writting data to disk')
            if count==0:
                with open(outFilePath,'w') as outFile:
                    np.savetxt(outFile,X=tempData)
                outFile=open(outFilePath,'a')
            else:
                np.savetxt(outFile,X=tempData)
            pbar.set_description('computing pair distances')
        count+=1
        gc.collect()
        pbar.update(chunkSize)
    pbar.close()
    outFile.close()
    
    if returnData:
        distData=np.concatenate(distData,axis=0)
        return(distData)

#Utilities to compute pearson and linear mutual information correlation
#matrices... these seem to be slow compared to carma and g_corr but may
#be useful when such tools cannot be easily compiled (e.g. cloud computing applications)
def calc_Ci(X,crdAxis=1):
    return(
        np.mean(
            np.apply_along_axis(
                lambda x: np.array(np.matrix([x]).T*np.matrix([x])),
                arr=X,
                axis=crdAxis),axis=0))

def calc_Cij(Xi,Xj,crdAxis=1):
    return(
        np.mean(
            np.apply_along_axis(
                lambda x: np.array(np.matrix([x]).T*np.matrix([x])),
                arr=np.concatenate([Xi,Xj],axis=crdAxis),
                axis=crdAxis),axis=0))

def calc_Linear_Mutual_Information(Xi,Xj,
                                   Ci=None,Cj=None,
                                   #Cii=None,Cjj=None,
                                   crdAxis=1,
                                   verbose=False):
    #Ci,Cii,Cj, and Cjj can be input if they have been calculated in advance
    #This can save significant when calcuting linear MI over all pairs in a 
    #large number of coordinate sets since Ci,Cii,Cj, and Cjj can be computed
    #in a single loop over all coordinate sets instead of needing to recalculated
    #for each ij coordinate set pair
    if Ci is None:
        ci=calc_Ci(Xi,crdAxis)
    else:
        ci=Ci
    #if Cii is None:
    #    cii=calc_Cij(Xi,Xi,crdAxis)
    #else:
    #    cii=Cii
    if Cj is None:
        cj=calc_Ci(Xj,crdAxis)
    else:
        cj=Cj
    #if Cjj is None:
    #    cjj=calc_Cij(Xj,Xj,crdAxis)
    #else:
    #    cjj=Cjj
        
    cij=calc_Cij(Xi,Xj,crdAxis)
    #cji=calc_Cij(Xj,Xi,crdAxis)
    CijMat=cij
    #CijMat=np.matrix(
    #    np.concatenate(
    #        [np.concatenate([cii,cij],axis=1),
    #         np.concatenate([cij,cjj],axis=1)],
    #        axis=0))
    if verbose:
        for entryName,entry in [
            ['Ci',ci], #['Cii',cii],
            ['Cj',cj], #['Cjj',cjj],
            #['Cij',cij],['Cji',cji],
            ['CijMat',CijMat],
            ['det(Ci)',np.linalg.det(ci)],
            ['det(Cj)',np.linalg.det(cj)],
            ['det(CijMat)',np.linalg.det(CijMat)]
        ]:
            print(entryName)
            print(entry)
    return(.5*(np.log(np.linalg.det(ci))+np.log(np.linalg.det(cj))-np.log(np.abs(np.linalg.det(CijMat)))))

def calc_pear_corr(Xi,Xj,Rii=None,Rjj=None,crdAxis=1,verbose=False):
    if Rii is None:
        rii=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xi**2,axis=crdAxis))
    else:
        rii=Rii
    if Rjj is None:
        rjj=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xj**2,axis=crdAxis))
    else:
        rjj=Rjj
    rij=np.mean(np.apply_along_axis(lambda x: np.sum(x),arr=Xi*Xj,axis=crdAxis))
    if verbose:
        print("rii:",rii)
        print("rjj:",rjj)
        print("rij:",rij)
    return(rij/(np.sqrt(rii)*np.sqrt(rjj)))
    

def calc_Normalized_LinearMI(Xi,Xj,
                             Ci=None,Cj=None,
                             #Cii=None,Cjj=None,
                             Rii=None,Rjj=None,
                             crdAxis=1,verbose=False):
    Imi=calc_Linear_Mutual_Information(Xi,Xj,
                                       Ci,Cj,
                                       #Cii,Cjj,
                                       crdAxis,
                                       verbose)
    Rmi=(1-np.exp(-2*Imi/Xi.shape[crdAxis]))**(1/2)
    
    rij=calc_pear_corr(Xi,Xj,Rii,Rjj,crdAxis)
    Igauss=-Xi.shape[crdAxis]/2.*np.log(1-rij**2)
    Rgauss=(1-np.exp(-2*Igauss/Xi.shape[crdAxis]))**(1/2)
    if verbose:
        print('rij',rij)
        print('Igauss',Igauss)
        print('Rgauss',Rgauss)
        print('Imi',Imi)
        print('Rmi',Rmi)
    return(Rmi)

def bokeh_interactive_matrix_plot(mat,nonzero_only=True,
                                    plotWidth=640,plotHeight=640,
                                    colorMap=sns.color_palette("coolwarm", n_colors=256).as_hex(),
                                    vmax=None,vmin=None,
                                    xlabel=None,ylabel=None,title=None
                                 ):
    if nonzero_only:
        matInds=np.nonzero(mat)
        matVals=np.array(mat[matInds]).flatten()
    else:
        matInds=(np.arange(mat.shape[0]),np.arange(mat.shape[1]))
        matVals=mat[matInds]
    
    tempFrame=pd.DataFrame({
        'X':matInds[0],
        'Y':matInds[1],
        'Val':matVals
    })
    
    vMax=np.max(matVals) if vmax is None else vmax
    vMin=np.min(matVals) if vmin is None else vmin
        
    valueColName='Val'
    
    mapper=LinearColorMapper(palette=colorMap,
                         low=vMin,
                         high=vMax)

    heatMapColumnSource=ColumnDataSource(tempFrame)

    toolTips=[(colName,"@{}".format(colName)) for colName in tempFrame.columns]

    p=bk.plotting.figure(plot_width=plotWidth,plot_height=plotHeight,
                tooltips=toolTips)

    if not (xlabel is None):
        p.xaxis.axis_label=xlabel
    if not (ylabel is None):
        p.yaxis.axis_label=ylabel
    if not (title is None):
        p.title.text=title
    
    p.background_fill_color="black"
    p.background_fill_alpha=.75
    
    color_bar=ColorBar(color_mapper=mapper,location='top_left')

    p.rect(x='X',y='Y',
           width=1,height=1,
           source=heatMapColumnSource,
           fill_color={'field':valueColName,'transform':mapper},line_color=None)
    p.add_layout(color_bar)
    show(p)
    
#A relatively simple bokeh based heatmap tool to quickly plot a heat map
#using a given value column and coordinate column pair from a data table
def bokeh_dataTable_heatMap(plotData,Xcol,Ycol,dataCol,
                            width=640,height=640,
                            rectheight=1,rectwidth=1,
                            colorMap=sns.color_palette("coolwarm", n_colors=256).as_hex(),
                            title=None,
                            xlabel=None,
                            ylabel=None,
                            axisFontSize="14pt",
                            vmin=None,vmax=None):
    
    if vmin is None:
        eMin=plotData[dataCol].min()
    else:
        eMin=vmin
    if vmax is None:
        eMax=plotData[dataCol].max()
    else:
        eMax=vmax
        
    plotData['ColorWeight']=plotData[dataCol] #.map(lambda x: cNorm(x)) #.map(np.log10)
    
    p=bk.plotting.figure(
        plot_width=width,plot_height=height)
    
    if title:
        p.add_layout(Title(text=title),'above')
    if xlabel:
        p.xaxis.axis_label=xlabel
    if ylabel:
        p.yaxis.axis_label=ylabel
        
    p.xaxis.major_label_text_font_size=axisFontSize
    p.yaxis.major_label_text_font_size=axisFontSize
    
    src=ColumnDataSource(plotData)
    mapper=linear_cmap(field_name='ColorWeight',
                       palette=colorMap,
                       low=eMin,
                       high=eMax,
                      )
    
    cbMapper=LinearColorMapper(palette=colorMap,
                               low=eMin,
                               high=eMax)
    color_bar = ColorBar(color_mapper=cbMapper, ticker= BasicTicker(),
                     location=(0,0))

    p.add_layout(color_bar, 'right')

    p.rect(source=src,width=rectwidth,height=rectheight, x=Xcol,y=Ycol,
             fill_color=mapper,color=mapper)

    
    
    hover=HoverTool()
    hover.tooltips=[
        (colName,"@{"+"{colName}".format(colName=colName)+"}")
        for colName in src.data.keys()
    ]
    print(hover.tooltips)
    p.add_tools(hover)

    output_notebook()
    show(p)
    
#functions for drawing networks using NGLview
import nglview as nv
import pytraj as pt

def bokeh_tidtyDataTable_barChart(dataTable,xColumns,valColumn,
                                  fillColumnInd=0,lineColumnInd=0,
                                  columnFormatFuns=None,
                                  xAxisLabel=None,yAxisLabel=None,
                                  pWidth=500,pHeight=900):
    plotData=dataTable[np.concatenate([xColumns,[valColumn]])].copy()
    if columnFormatFuns:
        for iCol,xCol in enumerate(xColumns):
            plotData[xCol]=plotData[xCol].map(columnFormatFuns[iCol])
    else:
        for xCol in xColumns:
            plotData[xCol]=plotData[xCol].map(str)
            
    xVals=[tuple(x) for x in list(plotData[xColumns].set_index(xColumns[0]).to_records())]
    
    source=ColumnDataSource(data={
        'x':xVals,
        valColumn:plotData[valColumn]
    })
    
    ToolTips=[(colName,'@%s'%colName) for colName in \
              source.to_df().columns]
    
    p=figure(x_range=FactorRange(*xVals),plot_height=500,plot_width=900,tooltips=ToolTips)
    p.vbar(x='x',top=valColumn,width=.9,source=source,
           fill_color=factor_cmap('x',palette=Spectral3,
                                  factors=plotData[xColumns[fillColumnInd]].unique(),
                                  start=fillColumnInd,end=fillColumnInd+1),
           line_color=factor_cmap('x',palette=Accent3,
                                  factors=plotData[xColumns[lineColumnInd]].unique(),
                                  start=lineColumnInd,end=lineColumnInd+1))
    
    p.xaxis.major_label_orientation = 87.5*np.pi/180.
    p.xaxis.group_label_orientation = 85.*np.pi/180.
    p.xaxis.subgroup_label_orientation = 80.*np.pi/180.
    if xAxisLabel:
        p.xaxis.axis_label = xAxisLabel
    if yAxisLabel:
        p.yaxis.axis_label = yAxisLabel
    show(p)

def drawProtNetEdge(protStruc,resID1,resID2,ngViewOb,
                    frame=0,edgeColor=[.5,.5,.5],radius=1,
                    *shapeArgs,**shapeKwargs):
    crd1=pt.center_of_mass(protStruc,':%g@CA'%resID1)[frame]
    crd2=pt.center_of_mass(protStruc,':%g@CA'%resID2)[frame]
    
    
    resname1=protStruc.topology.residue(resID1-1).name
    resid1=protStruc.topology.residue(resID1-1).original_resid
    
    resname2=protStruc.topology.residue(resID2-1).name
    resid2=protStruc.topology.residue(resID2-1).original_resid
    edgeLabel='%s.%g-%s.%g (%g-%g)'%(
        resname1,resid1,resname2,resid2,
        resID1-1,resID2-1)
    
    return ngViewOb.shape.add_cylinder(
        list(crd1),list(crd2),edgeColor,radius,
        edgeLabel,
        *shapeArgs,**shapeKwargs)

def getCorrNetEdgeColors(valMat,maskInds=None,
                         cmap=matplotlib.cm.get_cmap('viridis',1024),
                         cNorm=None,
                         baseColor=[0.,0.,0.]):
    if maskInds is None:
        nzInds=np.nonzero(np.isfinite(valMat))
    else:
        nzInds=maskInds
    
    if cNorm is None:
        vMin=np.min(valMat[nzInds])
        vMax=np.max(valMat[nzInds])
        norm=matplotlib.colors.Normalize(vmin=vMin,vmax=vMax)
    else:
        norm=cNorm
    
    tempCarray=cmap(norm(valMat[nzInds]).data)
    edgeColors=np.zeros((valMat.shape[0],valMat.shape[1],3))
    for ii in np.arange(3):
        edgeColors[:,:,ii]=baseColor[ii]
        edgeColors[:,:,ii][nzInds]=tempCarray[:,ii]
    return(edgeColors)

def drawEdgeColorBar(valMat,maskInds=None,
                     cmap=matplotlib.cm.get_cmap('viridis',1024),
                     cNorm=None,
                     ax=None,
                     barSize=[16,4],
                     label='Edge Color Scale',
                     label_size=10,
                     nTicks=5,
                     orientation='horizontal',
                     verbose=False):
    if maskInds is None:
        nzInds=np.nonzero(np.isfinite(valMat))
    else:
        nzInds=maskInds
    tempCmap=cmap
    tempCarray=tempCmap(valMat[nzInds])
    if verbose:
        print('vMin:',vMin,'vMax:',vMax)
    
    if cNorm is None:
        vMin=np.min(valMat[nzInds])
        vMax=np.max(valMat[nzInds])
        norm=matplotlib.colors.Normalize(vmin=vMin,vmax=vMax)
    else:
        norm=cNorm
    
    if ax is None:
        if orientation=='horizontal':
            fig = plt.figure(figsize=barSize)
            ax1 = fig.add_axes([0.3, 0.80, 0.4, 0.15])
        else:
            fig = plt.figure(figsize=np.flip(barSize))
            ax1 = fig.add_axes([0.3, 0.10, 0.15, 0.4])
    else:
        ax1=ax
    cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=tempCmap,
                                    norm=norm,
                                    orientation=orientation)
    cb1.set_label(label,size=label_size)
    cb1.ax.tick_params(labelsize=label_size)
    tick_locator = matplotlib.ticker.MaxNLocator(nbins=nTicks)
    cb1.locator = tick_locator
    cb1.update_ticks()
    plt.show()

def getCorrNetEdgeRadii(valMat,maskInds=None,
                        eMin=.0625,eMax=.75):
    if maskInds is None:
        nzInds=np.isfinite(valMat)
    else:
        nzInds=maskInds
        
    plotMat=valMat
    radiiMat=np.zeros(plotMat.shape)
    radiiMat[nzInds]=np.log(np.abs(np.array(plotMat)[nzInds]))
    rMin=np.min(radiiMat[nzInds])
    rMax=np.max(radiiMat[nzInds])
    rLinFun=lambda x: (x-rMin)/(rMax-rMin)
    rMap=lambda x:(eMax-eMin)*rLinFun(x)+eMin
    radiiMat[nzInds]=rMap(np.array(radiiMat[nzInds]))
    return(radiiMat)

def drawProtCorrMat(protStruc,corrMat,ngViewOb,
                    frame=0,colorsArray=None,radiiMat=None,
                    undirected=True):
    nzInds=np.nonzero(corrMat)
    edgeList=[]
    if not (radiiMat is None):
        radMat=radiiMat
    else:
        radMat=(np.abs(corrMat)>0)*.25
    if undirected:
        indArray=np.array([[ind[0],ind[1]] \
                            for ind in np.array([nzInds[0],nzInds[1]]).T \
                            if ind[0]<ind[1]])
    else:
        indArray=np.array([nzInds[0],nzInds[1]]).T
    for nzInd in indArray:
        if not (colorsArray is None):
            colorParm={'edgeColor':list(colorsArray[nzInd[0],nzInd[1],:])}
        else:
            colorParm={}
        edgeList.append(drawProtNetEdge(
            protStruc,nzInd[0]+1,nzInd[1]+1,
            ngViewOb,frame,radius=radMat[nzInd[0],nzInd[1]],
            **colorParm))
    return edgeList