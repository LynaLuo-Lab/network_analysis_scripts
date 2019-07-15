import numpy as np
import pandas as pd
import re
import sys
import os
import collections
import copy
import tqdm
from tqdm import tqdm
from tqdm import tqdm_notebook

noPytraj=False
try:
    import pytraj as pt
except ImportError:
    print "could not import pytraj"
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
            print "Reading file "+filePath+":",
        while line:
            lineTokens=line.split()
            lineEntries=np.extract(map(is_number,np.array(lineTokens)),np.array(lineTokens))
            entries[count:(count+len(lineEntries))]=lineEntries
            count=count+len(lineEntries)
            line=inputFile.readline()
            if verbose:
                if count % cMod == 0:
                    print str(np.floor(1000 * count / nEntries)/10)+"%",
        if verbose:
            print ""
    return collections.OrderedDict({"nRows":nRows,"nCols":nCols,"entries":entries})

def write_carma_matrix(filepath,nRows,nCols,dataEntries,
                       writeHeader=True,
                       header="((protein) and (not hydrogen)) and ((name CA) )",
                       verbose=False):
    if verbose:
        print "Writting carma matrix format file: "+filepath
    with open(filepath,'w') as outputFile:
        if writeHeader:
            if verbose:
                print 'Writting header'
            outputFile.write(header)
            outputFile.write('\n')
        if verbose:
            print 'Writting '+str(nRows)+' rows:',
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
                        print "\,    "
                        count=0
                    print "%4.1f%s "%((np.ceil((1000*iRow)/nRows)/10.0),"%"),
                    if nRows > 10000:
                        count = count+1
                else:
                    print ".",
        if verbose:
            print "\ndone"
            
            

def read_carma_matrix(filepath,has_header=True,returnHeader=False,verbose=False):
    #assumes all carma matrices are square matrices!
    if verbose:
        print "Reading carma matrix format file: "+filepath
    with open(filepath,'r') as inputFile:
        if(has_header):
            if verbose:
                print 'reading header'
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
                    print 'Reading '+str(nRows)+' rows: ',
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
                        print str(np.floor(1000*count/nRows)/10)+"% ",
                    else:
                        print ".",
                count=count+1
            line=inputFile.readline()
        if verbose:
            print ""
    if count > 0:
        if has_header & returnHeader:
            return collections.OrderedDict({"headerLine":headerLine,"nRows":nRows,"nCols":nRows,
                                            "entries":entries})
        else:
            return collections.OrderedDict({"nRows":nRows,"nCols":nRows,"entries":entries})
    else:
        print "Error! Data file appears empty!"

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
        print "Conversion complete"

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
        print "ERROR! data dictionary is missing needed key value pairs. !!Aborting!!"
        return
    if (writeHeader and useDictHeader):
        if not ('headerLine' in dataDict) :
            print "WARNING! headerLine was missing from data dictionary."
            print " -Defaulting to: "+header
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
        print "loading input trajectory"
    iterTraj=pt.iterload(inputTrajPath,top=inputTopPath)
    if verbose:
        print iterTraj
    commandList=map(lambda rid: "vector center :"+str(rid)+COM_mask,residList)
    if verbose:
        if len(commandList) >= 10:
            print "first 10 mapped commands: ",
            print commandList[0:10]
        else:
            print "mapped commands: ",
            print commandList
        print "running mapped commands"
    perResData=np.array(pt.compute(commandList,iterTraj,n_cores=threads).values())
    if resSelectMask=='':
        residString=map(lambda rid:str(rid),residList)
        resSelectMask=':'+','.join(residString)+resAtomMask
    if verbose:
        print "extracting residue selection subtrajectory"
    tempTraj=iterTraj[0:iterTraj.n_frames,resSelectMask]
    if verbose:
        print "subtrajectory info:"
        print tempTraj
    if not (tempTraj.shape[0]==perResData.shape[1] and \
            tempTraj.shape[1]==perResData.shape[0]):
        print "ERROR! mismatch between trajectory subselection coordinates and perResidue command data"
        print " trajectory.xyz.shape=",
        print trajectory.xyz.shape,
        print "; perResData.shape=",
        print perResData.shape[[1,0,2]]
        return
    if verbose:
        print "updating trajectory coordinates"
    for iDim in np.arange(3):
        tempTraj.xyz[:,:,iDim]=perResData[:,:,iDim].T
    print "done"
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
        print "Computing contact matrix:",
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
                print "%4.1f%s "%(np.ceil(1000*iRow/nResidues)/10,"%"),
                if lCount==10:
                    print "\n                         ",
                    lCount=0
                else:
                    lCount=lCount+1
        rgrid_i=residArray[mgrid_i[iRow,:].flatten()]
        rgrid_j=residArray[mgrid_j[iRow,:].flatten()]
        commandList=map(lambda rid,rjd: 'nativecontacts :'+\
                        str(rid)+' :'+str(rjd)+' mindist',rgrid_i,rgrid_j)
        tempDists=np.array(pt.compute(commandList,traj).values()[2:(3*len(rgrid_i)):3])
        if verbose and (verboseLevel>0):
            print commandList[0:3]
            print type(tempDists),
            print " ",
            print tempDists.shape,
            print tempDists
        matVals=map(lambda timeSeries:timeAggFun(timeSeries),
                    distSmoothFun(tempDists))
        tempMat[mgrid_i[iRow,:].flatten(),mgrid_j[iRow,:].flatten()]=matVals
        tempMat[iRow,iRow]=0
    print ""
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
        print "ERROR! Input matrices must by square!"
        eMat[:,:]=0
        return(eMat)
    if((Linv.shape[0]!=Amat.shape[0]) or
       (Linv.shape[1]!=Amat.shape[1])):
        print "ERROR! Input matrices must have the same shape!"
        eMat[:,:]=0
        return(eMat)
    if ((np.min(sources)<0) or
        (np.min(targets)<0) or
        (np.max(sources)>=Linv.shape[0]) or
        (np.max(targets)>=Linv.shape[1])):
        print "ERROR! invalid source or target index detected!"
        eMat[:,:]=0
        return(eMat)
    #get indices of edges... e.g. non-zero entries of Amat
    (Ei,Ej)=np.nonzero(Amat)
    if verbose:
        if verboseLevel > 2:
            print "Linv:",
            print Linv
            print "Amat:",
            print Amat
        print "computing betweenness for %g edges"%len(Ei)
        if verboseLevel > 0:
            print "Ei:",
            print Ei
            print "Ej:",
            print Ej
    if verbose and useProgressBar:
        Ebtw=np.array(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum(map(lambda src:
                        np.sum(map(lambda trg:
                            np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]),
                            targets)),
                        sources)),
                          pbarFun(Ei,file=sys.stdout),
                          Ej))/(len(sources)*len(targets))
    else:
        Ebtw=np.array(map(lambda i,j:
                    Amat[i,j]*\
                    np.sum(map(lambda src:
                        np.sum(map(lambda trg:
                            np.abs(Linv[i,src]+Linv[j,trg]-\
                                   Linv[i,trg]-Linv[j,src]),
                            targets)),
                        sources)),
                          Ei,
                          Ej))/(len(sources)*len(targets))
    if verbose:
        if verboseLevel > 0:
            print "Ebtw:",
            print Ebtw
            if verboseLevel > 1:
                print "(Ei,Ej):Ebtw;eMat"
    for iInd in np.arange(len(Ebtw)):
        eMat[Ei[iInd],Ej[iInd]]=Ebtw[iInd]
        if verbose :
            if verboseLevel > 1:
                print "(",
                print Ei[iInd],
                print ",",
                print Ej[iInd],
                print "):",
                print Ebtw[iInd],
                print ";",
                print eMat[Ei[iInd],Ej[iInd]]
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
        print "computing matrix Laplacian"
    Lmat=matLap(copy.deepcopy(mat))
    if verbose:
        print "extracting weighted adjacency matrix"
    Amat=matAdj(copy.deepcopy(mat))
    if verbose:
        print "computing moore-penrose inverse of matrix Laplacian"
    Linv=np.linalg.pinv(Lmat)
    if verbose:
        print "generating flow betweenness scores"
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
            print filepath+" does not exist."
        return False
    if not os.path.isfile(filepath):
        if verbose:
            print filepath+" is not a file."
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
            print filepath+' is a valid subopt file'
        return True
    else:
        if verbose:
            print filepath+" : ",
            if not foundPathStart:
                print "is missing paths section, ",
            if not foundPathCount:
                print "is missing path count"
        return False

def get_subopt_pathCount(filepath,verbose=False):
    if not validate_subopt_file(filepath,verbose=verbose):
        print "ERROR! "+filepath+" is not a valid subopt file."
        return -1
    else:
        with open(filepath,'r') as suboptFile:
            for line in suboptFile:
                if re.search('Number of paths is',line):
                    tokens=str.split(line)
                    pathCount=tokens[len(tokens)-1]
                    if verbose:
                        print 'path count = '+str(pathCount)
                    if not str.isdigit(pathCount):
                        break
                    return pathCount
        print "ERROR! Something went wrong, the file seemed valid but had an invalid path count line"
        return -1

def get_subopt_pathData(filepath,verbose=False):
    pathData=collections.OrderedDict()
    if not validate_subopt_file(filepath,verbose=verbose):
        print "ERROR! "+filepath+" is not a valid subopt file."
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
                        print 'path count = '+str(pathCount)
                    if not str.isdigit(pathCount):
                        foundPathCount=False
                        pathCount=-1
                    break
                if not foundPathStart:
                    if re.search('The final paths are',line):
                        foundPathStart=True
                else:
                    tokens=map(int,re.sub('[,)(]','',line).split())
                    tempPath=tokens[0:(len(tokens)-1)]
                    pathData['paths'].append(tempPath)
                    pathData['lengths'].append(tokens[len(tokens)-1])
        if not foundPathStart:
            print "ERROR! "+filepath+" seemed valid but path section was apparently absent"
            return pathData
        if not foundPathCount:
            print "Warning! final path count line was missing or ill-formed!"
            pathData['count']=len(pathData['paths'])
        else:
            if len(pathData['paths']) != int(pathCount):
                print "Warning! subopt file lists number of paths as",
                print str(pathCount),
                print "but",
                print len(pathData['paths']),
                print "paths were found."
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
        print 'running external subopt command: '+suboptCommand
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
        print 'running iterative dilation till '+str(minPathCount)+' paths are attained:'
        if not subVerbose:
            print '(%dilation,pathCount):',
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
        print "ERROR! subopt failed to generate a valid output file. Aborting"
        if returnSuboptData:
            if onlyFinalRun:
                return collections.OrderedDict()
            else :
                return [collections.OrderedDict()]
        else:
            return
    pathCount=get_subopt_pathCount(suboptDataFilePath)
    if verbose:
        print "(%g,%g)"%(float(percentDilation),float(pathCount)),
        if not subVerbose:
            print ",",
        else:
            print ""
    tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)
    optPathLength=np.min(map(float,tempData['lengths']))
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
            print "ERROR! subopt failed to generate a valid output file. Aborting"
            if returnSuboptData:
                return suboptoptData
            else:
                return
        pathCount=get_subopt_pathCount(suboptDataFilePath)
        if verbose:
            print "(%g,%g)"%(float(percentDilation),float(pathCount)),
        if not subVerbose:
            print ",",
        else:
            print ""
        tempData=get_subopt_pathData(suboptDataFilePath,subVerbose)   
        if returnSuboptData:
            if onlyFinalRun:
                suboptData=tempData
            else:
                suboptData.append(tempData)
    if verbose and not subVerbose:
        print ""
    if verbose:
        print "DONE!"
    if returnSuboptData:
        return suboptData
    
def get_subopt_dilations_data(suboptDir,basePattern,sep='.',
                              onlyMaxDilation=True,includeDilationValues=False,
                              includeFileNames=False,
                              verbose=False,verboseLevel=0):
    fileSearchPattern=str(sep).join([basePattern,'_'.join(['dilation','*']),'out'])
    if verbose:
        print 'file search pattern = '+fileSearchPattern
    searchPathPattern='/'.join([suboptDir,fileSearchPattern])
    if verbose:
        print 'searchPathPattern = '+searchPathPattern
    filePathList=glob.glob(searchPathPattern)
    fileNameList=map(lambda filepath: filepath.split('/')[-1],filePathList)
    if verbose and (verboseLevel > 0):
        print 'file name list: '+'\n'.join(fileNameList)
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
            print "Warning: "+suboptFileName+" was not a valid subopt file."
    if onlyMaxDilation:
        return suboptData
    else:
        return suboptDataSets

def get_index_of_nth_maxRank_element(valList,n,verbose=False):
    #this assumes that valList is already sorted in ascending order!
    u,v=np.unique(valList,return_inverse=True)
    maxRanks=(np.cumsum(np.bincount(v,minlength=u.size))-1)[v]
    if verbose:
        print 'value array:',
        print valList
        print 'rank array:',
        print maxRanks
    if np.max(maxRanks) < n:
        print "get nth maxRank: Warning! there are not enough elements"
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
        print 'max path index = '+str(maxPathIndex)
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
        map(lambda pathArray:'_'.join(map(str,pathArray)),pathData['paths']))

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
        print 'u:',
        print u
        print 'v:',
        print v
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
                  enumerate(netMatDict.keys()[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
    with tqdm(len(netMatDict)) as pbar:
        for matKey in netMatDict.keys():
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
                  enumerate(netMatDict.keys()[0].split(keySep))]
    else:
        keyNames=keyColNames
    if indexCols is None:
        indCols=keyNames
    else:
        indCols=indexCols
    with tqdm(len(netMatDict)) as pbar:
        for matKey in netMatDict.keys():
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