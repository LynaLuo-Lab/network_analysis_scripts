#require(rgl)
require(mgcv)
require(MASS)
require(ggplot2)
require(dplyr)
require(reshape2)
require(RColorBrewer)
require(colorRamps)
require(tidyr)
require(colorspace)
require(igraph)
require(SDMTools)

######################################DATA LOADING UTILITIES ######################################
#Individual data loading routines
#Ex use: 
# #Assumes you have the following in the directory path stored in variable 'dataFileDirectory'
# # 1) a carma matrix file 'contact.contact'
# # 2) a resinfo.cpptraj.dat file created using the 'resinfo' command in cpptraj
# # 3) a carma.average.dat file
# carmaData <- loard_carma_correlation_matrix(paste(dataFileDirectory,"contact.contact",sep="/"))
# resData <- load_system_info_data(dataFileDirectory)
# contactData <- apply_resData_to_carma_data(carmaData,resData)
load_carma_correlation_matrix <- function(datafilepath,has_carma_header=TRUE,
											indexNames=c("X","Y"),valuename="CORR",
											nnodes=NA) {
		#loads a single carma matrix as a dataframe
		#default frame will have colums 'X', 'Y', and 'CORR'
		if(has_carma_header==TRUE) { #skips first line, needed for most carma output files
			tempdatascan=scan(datafilepath,skip=1)
		} else {
			tempdatascan=scan(datafilepath)
		}
        if(is.na(nnodes)){ #guess dimensions from number of entries in matrix file
			datdims=c(sqrt(length(tempdatascan)),sqrt(length(tempdatascan)))
		} else { #explicitly set dimensions
			datadims=c(nnodes,nnodes)
		}
        tempdatarray<-array(tempdatascan,dim=datdims)
        tempdat <- tempdatarray %>% 
			melt(varames = c("X","Y"),value.name=valuename) %>% as.data.frame()
		tempdat[[indexNames[1]]]=tempdat[[colnames(tempdat)[1]]]
		tempdat[[indexNames[2]]]=tempdat[[colnames(tempdat)[2]]]
		tempdat<-tempdat %>% dplyr::select(4,5,3)
       	return(tempdat) 
}
load_system_info_data <- function(baseinfodir,
									resinfodir=baseinfodir,
									coordinfodir=baseinfodir,
									resinfofilename="resinfo.cpptraj.dat",
									coordfilename="carma.average.pdb") {
	#loads topology information from cpptraj atominof and resinfo
	#and structural information from carma pdb files into a data frame
	#defaults here assume all files are under the specified base directory (baseinfodir)
	resinfofile=paste(resinfodir,resinfofilename,sep="/")
	coordfile=paste(coordinfodir,coordfilename,sep="/")

	tempResDat=read.table(resinfofile,header=TRUE)

    tempResDat <- tempResDat %>% 
      mutate(Residue=paste(RESNAME,RESID,sep="_"))
    tempCoordDat <- read.table(coordfile,
                               col.names=c("Unused1","ResNum",
                                           "Unused2","RESNAME",
                                           "Chain","RESID",
                                           "rX","rY","rZ",
                                           "OCCUPANCY","BETA"),
                               fill = TRUE) %>% filter(Unused1 != "END")
    tempResDat <- left_join(tempResDat,
                            tempCoordDat %>% dplyr::select(ResNum,rX,rY,rZ),
                            by=c("ResNum"))

    resData <- tempResDat %>% 
                       dplyr::select(ResNum,RESID,RESNAME,Residue,rX,rY,rZ)
	return(resData)
}
apply_resData_to_carma_data <- function(carmadata,resdata,separator="_") {
	resInfoTypes=colnames(resdata)
	resInfoTypes=resInfoTypes[c(2:length(resInfoTypes))]
	carma_X_Name=colnames(carmadata[1])
	carma_Y_Name=colnames(carmadata[2])
	for(resInfoType in resInfoTypes) {
		X_InfoName=paste(carma_X_Name,resInfoType,sep=separator)
		Y_InfoName=paste(carma_Y_Name,resInfoType,sep=separator)
		carmadata[[X_InfoName]]=resdata[[resInfoType]][carmadata[[carma_X_Name]]]
		carmadata[[Y_InfoName]]=resdata[[resInfoType]][carmadata[[carma_Y_Name]]]
	}
	return(carmadata)
}
#Automated loading scripts to load data over an organized set of systems and windows
#SETUP information
#Directory structure diagram:
# -basedirectory-
# -> -system_directory- (name by "system_name")
# --> system_resinfo_file
# --> system_atominfo_file
# --> -window_directory- (name by "window_#")
# ---> contact_data_file
# ---> correlation_data_file
#
# corrleation and contact data file formats
# -- data files can be in either a table based format, or (now) original carma output format
# if the carma output format is used directly, set 'datafiles_as_arrays' and 
# 'contact_files_as_arrays' equal to "TRUE" respectively.
# if table based format is used:
# contact and correlation data must be three column (space separated) format
# first column should be labeled "X", second column "Y", and
# third column either "CORR" or "CONTACT" for correlation and contact data respectively
#
# system_list should be set as the list of systems (names of each 'system_directory')
# e.g. system_list=c("wt","R206H")
# window_list should be a list of window identifiers. Typically these should be integer
# values (other indexing may work... but do so at your own risk)
# e.g. window_list=c(1:10)
load_windowed_correlation_data <- function(basepath,baseinfodir=basepath,
                                           datafilename="network_varcov_table.dat",
                                           contactmapname="contact.contact.table.dat",
                                           resinfofilename="resinfo.cpptraj.dat",
                                           atominfofilename="atominfo.cpptraj.dat",
                                           coordfilename="carma.average.pdb",
                                           system_list,window_list,
										   skip_corr_header_rows=0,skip_contact_header_rows=1,
                                           datafiles_as_arrays=FALSE,
                                           contactfiles_as_arrays=FALSE,
                                           verbose=TRUE) {
  windowed_corr_data <- data.frame(SYSTEM=c(),WINDOW=c(),
                                   X=c(),Y=c(),
                                   X_RES=c(),Y_RES=c(),
                                   X_RESNAME=c(),Y_RESNAME=c(),
                                   X_Residue=c(),Y_Residue=c(),
                                   CORR=c())
  windowed_contact_map <- data.frame(SYSTEM=c(),WINDOW=c(),X=c(),Y=c(),
                                     X_RES=c(),Y_RES=c(),
                                     X_RESNAME=c(),Y_RESNAME=c(),
                                     X_Residue=c(),Y_Residue=c(),
                                     CONTACT=c())
  resData <- data.frame(SYSTEM=c(),
                        ResNum=c(),RESID=c(),RESNAME=c(),Residue=c(),
                        rX=c(),rY=c(),rZ=c())
  for(system in system_list) {
    if(verbose){print(paste("loading data for",system))
    print("-loading system structure data");flush.console()}
    sysinfodir=paste(baseinfodir,system,sep="/")
    tempResDat=read.table(paste(sysinfodir,resinfofilename,sep="/"),
                          header=TRUE)
    tempResDat <- tempResDat %>% mutate(SYSTEM=system) %>%
      dplyr::select(SYSTEM,ResNum,RESNAME,First,Last,RESID,MolNum)
    #read in atom information
    tempAtomDat=read.table(paste(sysinfodir,atominfofilename,sep="/"),
                           header=TRUE)
    tempAtomDat <- tempAtomDat %>% mutate(SYSTEM=system) %>%
      dplyr::select(SYSTEM,ATOMID,AtomNum,ATOMNAME,ResNum,RESNAME,MolNum,
                    ATOMTYPE,CHARGE,MASS,GBRADIUS,ELEMENT)
    #generate residue labels for plotting
    tempResDat <- tempResDat %>% 
      mutate(Residue=paste(RESNAME,RESID,sep="_"))
    tempCoordDat <- read.table(paste(basepath,system,"window_0",coordfilename,sep="/"),
                               col.names=c("Unused1","ResNum",
                                           "Unused2","RESNAME",
                                           "Chain","RESID",
                                           "rX","rY","rZ",
                                           "OCCUPANCY","BETA"),
                               fill = TRUE) %>% filter(Unused1 != "END")
    tempResDat <- left_join(tempResDat,
                            tempCoordDat %>% dplyr::select(ResNum,rX,rY,rZ),
                            by=c("ResNum"))
    resData <- rbind(resData,
                     tempResDat %>% 
                       dplyr::select(SYSTEM,ResNum,RESID,RESNAME,Residue,rX,rY,rZ))
    #add residue labels to atom information data
    tempAtomDat <- tempAtomDat %>% 
      mutate(RESID=tempResDat$RESID[ResNum],
             Residue=tempResDat$Residue[ResNum])
    #sanity check on ATOMID column in wavelet data
    resMax=length(tempResDat$RESID)
    if(verbose){print("-loading windowed correlation and contact data:");cat("-   ");flush.console()}
    for(window in window_list) {
      if(verbose){cat(paste(" ",window));if(window%%25 == 24){cat("/n-  ")};flush.console()}
      datafilepath=paste(basepath,"/",system,"/window_",window,"/",datafilename,sep="")
      #filter
      contactpath=paste(basepath,"/",system,"/window_",window,"/",contactmapname,sep="")
      if (datafiles_as_arrays) {
        #tempdatscan=scan(datafilepath,skip=skip_corr_header_rows)
        #datdims=c(sqrt(length(tempdatscan)),sqrt(length(tempdatscan)))
        #tempdatarray<-array(tempdatscan,dim=datdims)
        #tempdat <- tempdatarray %>% melt(varames = c("X","Y"),value.name="CORR")
        #rm(tempdatscan)
        #rm(datdims)
        #rm(tmpdatarray)
		tempdat <- load_carma_correlation_matrix(datafilepath,has_carma_header=(skip_corr_header_rows > 0),
											indexNames=c("X","Y"),valuename="CORR")
      } else {
        tempdat=read.table(datafilepath,header = TRUE)
      }
      if (contactfiles_as_arrays) {
        #tempcontactscan=scan(contactpath,skip=skip_contact_header_rows)
        #contactdims=c(sqrt(length(tempcontactscan)),sqrt(length(tempcontactscan)))
        #tempcontactarray<-array(tempcontactscan,dim=contactdims)
        #tempcontact<-tempcontactarray %>% melt(varnames=c("X","Y"),value.name="CONTACT")
        #rm(tempcontactscan)
        #rm(contactdims)
        #rm(tempcontactarray)
		tempcontact <- load_carma_correlation_matrix(datafilepath,has_carma_header=(skip_contact_header_rows > 0),
											indexNames=c("X","Y"),valuename="CONTACT")
      } else {
        tempcontact=read.table(contactpath,header = TRUE)
      }
      tempcontact <- tempcontact %>%
        mutate(SYSTEM=system,WINDOW=window,
               X_RES=tempResDat$RESID[X],
               Y_RES=tempResDat$RESID[Y],
               X_RESNAME=tempResDat$RESNAME[X],
               Y_RESNAME=tempResDat$RESNAME[Y],
               X_Residue=tempResDat$Residue[X],
               Y_Residue=tempResDat$Residue[Y]
        ) %>%
        dplyr::select(SYSTEM,WINDOW,
                      X,Y,X_RES,Y_RES,
                      X_RESNAME,Y_RESNAME,X_Residue,Y_Residue,
                      CONTACT)
      #tempdat <- left_join(tempdat,tempcontact,by=c("X","Y"))
      tempdat <- tempdat %>% 
        #filter(CONTACT==1) %>% #filters out correlations not in contact map
        #dplyr::select(X,Y,CORR) %>% #remove contact map column
        filter(X>0,X<resMax,Y>0,Y<resMax) %>%
        mutate(X_RES=tempResDat$RESID[X],
               Y_RES=tempResDat$RESID[Y],
               X_RESNAME=tempResDat$RESNAME[X],
               Y_RESNAME=tempResDat$RESNAME[Y],
               X_Residue=tempResDat$Residue[X],
               Y_Residue=tempResDat$Residue[Y],
               SYSTEM=system,WINDOW=window) %>%
        dplyr::select(SYSTEM,WINDOW,
                      X,Y,
                      X_RES,Y_RES,
                      X_RESNAME,Y_RESNAME,
                      X_Residue,Y_Residue,
                      CORR)
      windowed_corr_data=rbind(
        windowed_corr_data,
        tempdat
      )
      windowed_contact_map <- rbind(windowed_contact_map,tempcontact)
    }
    if(verbose){cat("\n");flush.console()}
  }
  if(verbose){print("Done!");flush.console()}
  return(list(windowed_corr_data=windowed_corr_data,
              windowed_contact_data=windowed_contact_map,
              res_data=resData))
}

#merge windowed correlation data and windowed contact data
gen_contact_corr_summary <- function(windowed_corr_data,windowed_contact_data,
                                     verbose=FALSE,vlevel=1,
									 calc_pvals=FALSE,use_abs_corr=FALSE) {
  if (verbose) {print("melting windowed_corr_data");flush.console()}
  melted_windowed_corr <- melt(windowed_corr_data,
                               id.vars=c("SYSTEM","WINDOW","X","Y",
                                         "X_RES","Y_RES","X_RESNAME","Y_RESNAME",
                                         "X_Residue","Y_Residue"),
                               measure.vars=c("CORR"),
                               variable.name="measurement",
                               value.name="value")
  if (verbose) {print("melting windowed_contact_data");flush.console()}
  melted_windowed_contact <- melt(windowed_contact_data,
                                  id.vars=c("SYSTEM","WINDOW","X","Y",
                                            "X_RES","Y_RES","X_RESNAME","Y_RESNAME",
                                            "X_Residue","Y_Residue"),
                                  measure.vars=c("CONTACT"),
                                  variable.name="measurement",
                                  value.name="value")
  if (verbose) {print("joining melted data");flush.console()}
  melted_windowed_contact_corr <- rbind(melted_windowed_corr,melted_windowed_contact)
  if(verbose && vlevel>1){melted_windowed_contact %>% glimpse;flush.console()}
  rm(melted_windowed_corr)
  rm(melted_windowed_contact)
  if (verbose) {print("re-spreading melted data");flush.console()}
  windowed_contact_corr <- melted_windowed_contact_corr %>%
    tidyr::spread(key=measurement,value=value,fill=0.0)
  rm(melted_windowed_contact_corr)
  if (verbose) {print("computing windowed_data summary");flush.console()}
  if (use_abs_corr) {
    windowed_contact_corr <- windowed_contact_corr %>%
	  mutate(CORR = abs(CORR))
  }
  if (!calc_pvals) {
    contact_corr_summary <- windowed_contact_corr %>%
      group_by(SYSTEM,X,Y,X_RES,Y_RES,X_RESNAME,Y_RESNAME,X_Residue,Y_Residue) %>%
      summarise(CORR_CONTACT_sum=sum(CORR*CONTACT),
                CONTACT_sum=sum(CONTACT),
				CONTACT_stdev=sd(CONTACT),
                CONTACT=mean(CONTACT)) %>%
      filter(CONTACT>0) %>% mutate(CORR=CORR_CONTACT_sum/CONTACT_sum) %>%
      dplyr::select(SYSTEM,X,Y,X_RES,Y_RES,X_RESNAME,Y_RESNAME,X_Residue,Y_Residue,
                    CONTACT,CORR,CONTACT_stdev) %>%
      as.data.frame()
  } else {
    contact_corr_summary <- windowed_contact_corr %>%
      group_by(SYSTEM,X,Y,X_RES,Y_RES,X_RESNAME,Y_RESNAME,X_Residue,Y_Residue) %>%
      summarise(CORR_CONTACT_sum=sum(CORR*CONTACT),
                CONTACT_sum=sum(CONTACT),
				CORR_pValue=t.test(CONTACT*CORR)$p.value,
				CONTACT_stdev=sd(CONTACT),
                CONTACT=mean(CONTACT)
			   ) %>%
      filter(CONTACT>0) %>% mutate(CORR=CORR_CONTACT_sum/CONTACT_sum) %>%
      dplyr::select(SYSTEM,X,Y,X_RES,Y_RES,X_RESNAME,Y_RESNAME,X_Residue,Y_Residue,
                    CONTACT,CORR,
					CONTACT_stdev,
					CORR_pValue) %>%
      as.data.frame()
  }
  return(contact_corr_summary)
}

#Loads tables of flow-betweenness data and corresponding subOptimal paths (if loadSubOptData=TRUE)
#for given sets of sources and targets.
#For the subOptimal paths, only the edge names are saved. These must correspond exactly with
#the edge names in the assosciated network data table.
#The network data tables must have the columns:
#( name from to weight CORR CONTACT rDist btw )
#These should be included as column names in the header
#the columns "name" "from" and "to" should be the edges of a correlation network graph
#with CORR describing the correlation, and btw the betweenness scores of each edge
#for the given source-target set pair
#The subOpt data table should follow the same structure, but only the edge names are truly
#needed. The other data is discarded.
#The return value will be a list of lists (3 layers deep) with the outermost list indexed by
#system and the next sublist will contain the vertex data for the system in an entry
#labeled vertexData along with lists of edge and subopt data indexed by source-target pair.
#The innermost list will contain the entries
#"network" and (if loadSubOptData=TRUE) "subOptEdges"
#which contain the correlation network and array of subOptimal path edges respectively
load_corr_flow_btw_tables <- function(basepath,baseinfodir=basepath,
                                      networkDataExt="NetworkTable.dat",
                                      subOptDataExt="SubOptTable.dat",
                                      resinfofilename="resinfo.cpptraj.dat",
                                      atominfofilename="atominfo.cpptraj.dat",
                                      coordfilename="carma.average.pdb",
                                      system_list,src_list,trg_list,loadSubOptData=TRUE,
                                      verbose=TRUE) {
  dataList=list()
  resData <- data.frame(SYSTEM=c(),
                        ResNum=c(),RESID=c(),RESNAME=c(),Residue=c(),
                        rX=c(),rY=c(),rZ=c())
  for(system in system_list) {
    if(verbose){print(paste("loading data for",system))
      print("-loading system structure data");flush.console()}
    sysinfodir=paste(baseinfodir,system,sep="/")
    tempResDat=read.table(paste(sysinfodir,resinfofilename,sep="/"),
                          header=TRUE)
    tempResDat <- tempResDat %>% mutate(SYSTEM=system) %>%
      dplyr::select(SYSTEM,ResNum,RESNAME,First,Last,RESID,MolNum)
    #read in atom information
    tempAtomDat=read.table(paste(sysinfodir,atominfofilename,sep="/"),
                           header=TRUE)
    tempAtomDat <- tempAtomDat %>% mutate(SYSTEM=system) %>%
      dplyr::select(SYSTEM,ATOMID,AtomNum,ATOMNAME,ResNum,RESNAME,MolNum,
                    ATOMTYPE,CHARGE,MASS,GBRADIUS,ELEMENT)
    #generate residue labels for plotting
    tempResDat <- tempResDat %>% 
      mutate(Residue=paste(RESNAME,RESID,sep="_"))
    tempCoordDat <- read.table(paste(basepath,system,"window_0",coordfilename,sep="/"),
                               col.names=c("Unused1","ResNum",
                                           "Unused2","RESNAME",
                                           "Chain","RESID",
                                           "rX","rY","rZ",
                                           "OCCUPANCY","BETA"),
                               fill = TRUE) %>% filter(Unused1 != "END")
    tempResDat <- left_join(tempResDat,
                            tempCoordDat %>% dplyr::select(ResNum,rX,rY,rZ),
                            by=c("ResNum"))
    resData <- rbind(resData,
                     tempResDat %>% 
                       dplyr::select(SYSTEM,ResNum,RESID,RESNAME,Residue,rX,rY,rZ))
    #add residue labels to atom information data
    tempAtomDat <- tempAtomDat %>% 
      mutate(RESID=tempResDat$RESID[ResNum],
             Residue=tempResDat$Residue[ResNum])
    #sanity check on ATOMID column in wavelet data
    resMax=length(tempResDat$RESID)
    if(verbose){print("-loading network correlation and betweenness data:");cat("-   ");flush.console()}
    for(src in src_list) {
      for(trg in trg_list) {
        dataPrefix=paste(system,src,trg,sep="_")
        if(verbose){
          print(paste("--loading data for",
                      dataPrefix,"correlation flow network and sub-Optimal paths"))
          flush.console()
        }
        networkFilePath=paste(basepath,"/",system,"/",dataPrefix,".",networkDataExt,sep="")
        subOptFilePath=paste(basepath,"/",system,"/",dataPrefix,".",subOptDataExt,sep="")
        tempNetwork <- read.table(networkFilePath,header=TRUE,stringsAsFactors=TRUE)
        tempNetwork <- tempNetwork %>%
          mutate(X=from,Y=to)
        tempCols=colnames(tempNetwork)
        tempNetwork <- tempNetwork %>%
          dplyr::select(X,Y,name,weight,CORR,CONTACT,rDIST,btw) %>%
          mutate(SYSTEM=system,SOURCE=src,TARGET=trg,NetName=dataPrefix,
                 X_RES=tempResDat$RESID[X],
                 Y_RES=tempResDat$RESID[Y],
                 X_RESNAME=tempResDat$RESNAME[X],
                 Y_RESNAME=tempResDat$RESNAME[Y],
                 X_Residue=tempResDat$Residue[X],
                 Y_Residue=tempResDat$Residue[Y]
          ) %>%
          dplyr::select(SYSTEM,SOURCE,TARGET,NetName,
                        X,Y,name,X_RES,Y_RES,
                        X_RESNAME,Y_RESNAME,X_Residue,Y_Residue,
                        CONTACT,CORR,weight,rDIST,btw)
        
        if (loadSubOptData) {
          subOptData <- read.table(subOptFilePath,header=TRUE,stringsAsFactors=TRUE)
          dataList[[system]][[paste(src,trg,sep="_")]]=list(
            network=tempNetwork,subOptEdges=subOptData$name
          )
          dataList[[system]][["vertexData"]]=tempResDat
        } else {
          dataList[[system]][[paste(src,trg,sep="_")]]=list(network=tempNetwork)
          dataList[[system]][["vertexData"]]=tempResDat
        }
        
      }
    }
  }
  if(verbose){print("Done!");flush.console()}
  return(dataList)
}

#convenience functions for use in forcing matrix symmetry
mergeWhichNonZeroAbsMin=function(x,y){
	#this function has a somewhat more complicated structure
	#since it will most likely be used to write weight matrices
	#which use 0 to indicate the absence of an edge, but otherwise
	#use the inverse log of a corresponding value.
	#E.g. low values result in very high weights and vice-versa
	#Thus if both entires are zero should be used if both
	#entries are non-finite. This is useful since Inf, or NA
	#may be commone if -log() is applied to an entry where the
	#the value was zero. In that case, we probably want to turn
	#the edge off.
	if((is.finite(x) & is.finite(y))){
		if((abs(x)>0) & (abs(y)>0)) {
			return(c(x,y)[which.min(abs(c(x,y)))])
		} else {
			if(abs(x)>0) {
				return(x)
			} else if (abs(y) > 0) {
				return(y)
			} else {
				return(0)
			}
		}
	} else {
		if(is.finite(x)) {
			return(x)
		} else if(is.finite(y)) {
			return(y)
		} else {
			return(0)	
		}
	}
}
mergeWhichAbsMax=function(x,y){
    return(c(x,y)[which.max(abs(c(x,y)))])
}
mergeWhichAbsMin=function(x,y){
    return(c(x,y)[which.min(abs(c(x,y)))])
}
mergeAbsMax=function(x,y){
    return(max(abs(c(x,y))))
}
mergeAbsMin=function(x,y){
    return(min(abs(c(x,y))))
}
mergeMax=function(x,y){
    return(max(c(x,y)))
}
mergeMin=function(x,y){
    return(min(c(x,y)))
}
mergeMean=function(x,y){
    return(mean(c(x,y)))
}
mergeAbsMean=function(x,y){
	return(mean(abs(c(x,y))))
}

makeSymmetric=function(mat,mergeFun=mergeWhichAbsMax){
    tmat=t(mat)
    outMat=(mat+tmat)*0
    dimlim=min(dim(mat))
    for(ii in c(1:(dimlim-1))){
        for(jj in c(ii:dimlim)){
            val=mergeFun(mat[ii,jj],mat[jj,ii])
            outMat[ii,jj]=val
            outMat[jj,ii]=val
        }
    }
    return(outMat)
}

symmetricSparseMatrixFromLists <- function(fromList,toList,valList,n,
									 mergeFunction=mergeWhichAbsMax) {
	tempMat=sparseMatrix(i=fromList,j=toList,x=valList,dims=c(n,n))
	outMat=tempMat+t(tempMat)
	for(iInd in 1:length(fromList)){
		ii=fromList[iInd]
		jj=toList[iInd]
		val=mergeFunction(tempMat[ii,jj],tempMat[jj,ii])
		outMat[ii,jj]=val
		outMat[jj,ii]=val
	}
	return(outMat)
}

write_networkView_matrix <- function(outputFileName="data.dat",
                                     header="((protein) and (not hydrogen)) and ((name CA) )",
                                     fromList,toList,valList,
                                     nodeCount=length(unique(c(fromList,toList))),
                                     makeSymmetric=FALSE,
									 mergeFunction=mergeWhichAbsMax,
									 verbose=FALSE){
  if(makeSymmetric) {
	tempMat=symmetricSparseMatrixFromLists(fromList,toList,valList,nodeCount,mergeFunction)
  } else {
	tempMat=sparseMatrix(i=fromList,j=toList,x=valList,dims=c(nodeCount,nodeCount))
  }
    if(verbose) {
	toColumn=toList
	fromColumn=fromList
	valColumn=valList
  } 

  if(verbose) {
	message("-NetworkView Matrix Data Summary-")
	message("")
	message(" Input: Min 1stQt. Median Mean 3rdQt. Max")
	message(" toColumn: ",(toColumn) %>% summary %>% paste(collapse=" "))
	message(" fromColumn: ",(fromColumn) %>% summary %>% paste(collapse=" "))
	message(" valColumn: ",(valColumn) %>% summary %>% paste(collapse=" "))
	message("")
	message(" Output: Min 1stQt. Median Mean 3rdQt. Max")
	message(" toColumn: ",(tempMat %>% summary)[["i"]] %>% summary %>% paste(collapse=" "))
	message(" fromColumn: ",(tempMat %>% summary)[["j"]] %>% summary %>% paste(collapse=" "))
	message(" valColumn: ",(tempMat %>% summary)[["x"]] %>% summary %>% paste(collapse=" "))
	message("")
	message("-",rep("#",31)%>%paste(collapse=""),"-")
  }

  tempMat=tempMat%>%as.matrix()
  outFile=file(outputFileName,"w")
  cat(header,file=outFile)
  cat("\n",file=outFile)
  for (iRow in c(1:nrow(tempMat))){
    for (iCol in c(1:ncol(tempMat))){
      if(abs(tempMat[iRow,iCol]) >= 1E-7) {
        cat(sprintf("%9.5f ",tempMat[iRow,iCol]),file=outFile)
      } else {
        cat("0 ",file=outFile)
      }
    }
    # cat(tempMat[iRow,],file=outFile)
    cat("\n",file=outFile)
  }
  close(outFile)
}
write_networkDataTable_to_networkView_matrix <- function(outputFileName="data.dat",
                                                         header="((protein) and (not hydrogen)) and ((name CA) )",
                                                         networkData,
                                                         toColumn="X",fromColumn="Y",
                                                         valColumn="CORR",
                                                         nodeCount=0,
                                                         makeSymmetric=FALSE,mergeFunction=mergeWhichAbsMax,
														 indexOffset=0,verbose=FALSE) {
  if(verbose) {
	message("-NetworkView Matrix Data Table Summary-")
	message(" Column: Min 1stQt. Median Mean 3rdQt. Max")
	message(" toColumn: ",(networkData[[toColumn]]-indexOffset) %>% summary %>% paste(collapse=" "))
	message(" fromColumn: ",(networkData[[fromColumn]]-indexOffset) %>% summary %>% paste(collapse=" "))
	message(" valColumn: ",(networkData[[valColumn]]) %>% summary %>% paste(collapse=" "))
	message("-",rep("#",37)%>%paste(collapse=""),"-")
  }
  fromInds=networkData[[toColumn]]-indexOffset
  toInds=networkData[[fromColumn]]-indexOffset
  xVals=networkData[[valColumn]]
  if (nodeCount==0) {
    nodeCount=unique(c(fromInds,toInds))
  }
  write_networkView_matrix(outputFileName,header,fromInds,toInds,xVals,nodeCount=nodeCount,
	makeSymmetric=makeSymmetric,mergeFunction=mergeFunction,verbose=verbose)
}
#routines for loading suboptimal path data generated by the VMD NetworkView plugin 'subopt' command
read_networkView_suboptFile_pathCount <- function(suboptFilePath,verbose=FALSE,verboseLevel=0) {
	tempFile=file(suboptFilePath,"r")
	tempLine=readLines(tempFile,n=1)
	if(length(tempLine)==0) {
		print("Error: Data file appears to be empty!")
		close(tempFile)
		return(0)
	} else {
		linenum=1
		if(verbose && verboseLevel>1){print(paste("- line number",linenum," = ",tempLine))}
		startLineFound=tempLine %>% grepl("Number of paths is",.)
		while(!startLineFound) {
			tempLine=readLines(tempFile,n=1)
			linenum=linenum+1
			if(length(tempLine)==0) {
				if(verbose && verboseLevel>0){print(paste("- line number",linenum,"was empty"))}
				break
			} else {
				if(verbose && verboseLevel>1){print(paste("- line number",linenum," = ",tempLine))}
			}
			startLineFound=tempLine%>%grepl("Number of paths is",.)
		}
		if(!startLineFound){
			print("Error: Data file ended prematurely or has wrong format!")
			close(tempFile)
			return(0)
		} else {
			lineTokens <- (tempLine%>%strsplit(" "))[[1]] %>% as.character
			pathCount=lineTokens[length(lineTokens)]%>%as.numeric
			close(tempFile)
			return(pathCount)
		}
	}
}
read_networkView_suboptFile_nodeCounts <- function(suboptFilePath,resCount,verbose=FALSE,verboseLevel=1) {
  tempFile=file(suboptFilePath,"r")
  tempLine=readLines(tempFile,n=1)
  tempMat=c(rep(0,resCount))
  if(length(tempLine)==0){
    print("Error: Data File appears to be empty!")
    close(tempFile)
    return(tempMat)
  } else {
    startLineFound=tempLine %>% grepl("The final paths are",.)
    while(!startLineFound) {
      tempLine=readLines(tempFile,n=1)
      if (length(tempLine)==0) {
        break
      }
      startLineFound=tempLine %>% grepl("The final paths",.)
    }
    if (!startLineFound) {
      print("Error: Data File ended prematurely or has wrong format!")
      close(tempFile)
      return(tempMat)
    } else {
      #print(tempLine)
      tempLine=readLines(tempFile,n=1)
      if(length(tempLine)==0) {
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      }
      endLineFound=tempLine %>% grepl("Number of paths is",.)
      if (endLineFound){
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      } else {
        count=1
        #tempMat=sparseMatrix(i=c(1),j=c(1),x=c(0),dims=c(resCount,resCount))
        while(!endLineFound) {
          if(verbose){
            print(paste("**parsing subopt path",count))
            if(verboseLevel>1) {
              print(paste("***line = ",tempLine))
            }
          }
          tempVals=(tempLine%>%strsplit("[ ,()]+"))[[1]] %>% as.numeric()
          tempVals=tempVals+1
          if(verbose & verboseLevel > 2) {
            print(paste("****line tokens = [",
                        paste(tempVals,collapse = ", "),
                        "]"))
          }
          tempMat[tempVals[1:(length(tempVals)-1)]]=tempMat[tempVals[1:(length(tempVals)-1)]]+1 
          tempLine=readLines(tempFile,n=1)
          if (length(tempLine)==0) {
            break
          }
          endLineFound=tempLine %>% grepl("Number of paths is",.)
          count=count+1
        }
      }
      if (!endLineFound) {
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      } else {
        close(tempFile)
        return(tempMat)
      }
    }
  }

}
read_networkView_suboptFile_to_Matrix <- function(suboptFilePath,resCount,verbose=FALSE,verboseLevel=1) {
  tempFile=file(suboptFilePath,"r")
  tempLine=readLines(tempFile,n=1)
  tempMat=sparseMatrix(i=c(1),j=c(1),x=c(0),dims=c(resCount,resCount))
  if(length(tempLine)==0){
    print("Error: Data File appears to be empty!")
    close(tempFile)
    return(tempMat)
  } else {
    startLineFound=tempLine %>% grepl("The final paths are",.)
    while(!startLineFound) {
      tempLine=readLines(tempFile,n=1)
      if (length(tempLine)==0) {
        break
      }
      startLineFound=tempLine %>% grepl("The final paths",.)
    }
    if (!startLineFound) {
      print("Error: Data File ended prematurely or has wrong format!")
      close(tempFile)
      return(tempMat)
    } else {
      #print(tempLine)
      tempLine=readLines(tempFile,n=1)
      if(length(tempLine)==0) {
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      }
      endLineFound=tempLine %>% grepl("Number of paths is",.)
      if (endLineFound){
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      } else {
        count=1
        #tempMat=sparseMatrix(i=c(1),j=c(1),x=c(0),dims=c(resCount,resCount))
        while(!endLineFound) {
          if(verbose){
            print(paste("**parsing subopt path",count))
            if(verboseLevel>1) {
              print(paste("***line = ",tempLine))
            }
          }
          tempVals=(tempLine%>%strsplit("[ ,()]+"))[[1]] %>% as.numeric()
          tempVals=tempVals+1
          if(verbose & verboseLevel > 2) {
            print(paste("****line tokens = [",
                        paste(tempVals,collapse = ", "),
                        "]"))
          }
          nVals=length(tempVals)
          nNodes=nVals-1
          pathLength=tempVals[nVals]
          tempPathMat=sparseMatrix(i=tempVals[1:(nNodes-1)]%>%as.numeric(),
                                   j=tempVals[2:nNodes]%>%as.numeric(),
                                   x=rep(1,nNodes-1),
                                   dims=c(resCount,resCount))
          for(iInd in c(1:(nNodes-1))) { #ensure the matrix is symmetric
            ii=tempVals[iInd]
            jj=tempVals[iInd+1]
            tempPathMat[jj,ii]=1
            tempPathMat[ii,jj]=1
          }
          #tempPathMat=tempPathMat+t(tempPathMat)
          tempMat=tempMat+tempPathMat
          tempLine=readLines(tempFile,n=1)
          if (length(tempLine)==0) {
            break
          }
          endLineFound=tempLine %>% grepl("Number of paths is",.)
          count=count+1
        }
      }
      if (!endLineFound) {
        print("Error: Data File ended prematurely or has wrong format!")
        close(tempFile)
        return(tempMat)
      } else {
        close(tempFile)
        return(tempMat)
      }
    }
  }
}
read_networkView_suboptFile_nodeCounts_to_DataFrame <- function(suboptFilePath,resCount) {
  tempMat=read_networkView_suboptFile_nodeCounts(suboptFilePath,resCount)
  if(!is.null(tempMat) && 
     !is.null(length(tempMat) == 0)) {
    tempFrame=data.frame(ResNum=c(1:resCount),COUNT=tempMat)
  } else {
    tempFrame=data.frame(ResNum=c(),COUNT=c())
  }
  return(tempFrame)
}

read_networkView_suboptFile_to_DataFrame <- function(suboptFilePath,resCount) {
  tempMat=read_networkView_suboptFile_to_Matrix(suboptFilePath,resCount)
  if(!is.null(tempMat) && 
     !is.null(dim(tempMat)) && 
     dim(tempMat)[1]==resCount &&
     dim(tempMat)[2]==resCount){
    tempFrame=tempMat%>%as.matrix%>%melt%>%as.data.frame()%>%
      rename(X=Var1,Y=Var2,COUNT=value)%>%filter(COUNT>0)
  } else {
    tempFrame=data.frame(X=c(),Y=c(),COUNT=c())
  }
  return(tempFrame)
}

load_windowed_networkView_suboptData <- function(suboptDataDir,suboptFileExt,
                                        systemNames,windows,weightTypes,resCountsList,
                                        sourceSetNamesList,sourceSetNodeIDsList,
                                        targetSetNamesList,targetSetNodeIDsList,
                                        verbose=FALSE,verboseLevel=0,ordered_nodes=FALSE) {
  #reads sets of subopt files created by networkView into a joint data frame
  #
  #all subopt data files should be stored in suboptDataDir and named as:
  #	systemName_sourceSetName_targetSetName.window_windowID.weightType.sourceNode_targetNode.suboptFileExt
  # e.g. wt_FKBP_res375.window_1.CORRweight.203_375.out.out
  #
  #Argument specifications:
  ### systemNames: 			a vector of strings representing the names of each system to load
  ###							subopt network data for
  ###
  ### windows:				a list / array of window numbers or names
  ### weightTypes: 			a vector of strings representing the names of each weighting used to
  ###				  			calculate different subopt path networks
  ###
  ### resCountsList:			list containing the number of residues for each system
  ###							i.e. accessed as: resCountsList[[systemName]]
  ###
  ### sourceSetNamesList:		a 'list()' of vectors of source set names organized by system name
  ###							i.e. accessed as: sourceSetNamesList[[systemName]]
  ###
  ### sourceSetNodeIDsList:	a list of lists 'list(list())' of vectors of source node indices
  ###							organized by system name then source set name.
  ###							i.e. accesssed as: sourceSetNodeIDsList[[systemName]][[sourceSetName]]
  ### targetSetNamesList:		same as sourceSetNamesList, but for target sets
  ### 
  ### targetSetNodeIDsList:	same as sourceSetNodeIDsList, but of target sets
  suboptDataFrame=data.frame(SYSTEM=c(),WINDOW=c(),WEIGHT_TYPE=c(),
                             SOURCE_SET=c(),TARGET_SET=c(),
                             SOURCE_NODE=c(),TARGET_NODE=c(),
							 PATHS_COUNT=c(),
                             X=c(),Y=c(),COUNT=c())
  for(systemName in systemNames){
    if(verbose){print(paste("working on system",systemName))}
    sourceSetNames=sourceSetNamesList[[systemName]]
    targetSetNames=targetSetNamesList[[systemName]]
    resCount=resCountsList[[systemName]]
	for(window in windows) {
		if(verbose){print(paste("-working on window",window))}
    	for(sourceSetName in sourceSetNames){
    	  sourceSetNodeIDs=sourceSetNodeIDsList[[systemName]][[sourceSetName]]
    	  for(targetSetName in targetSetNames){
    	    if(verbose && verboseLevel>0){print(paste("--working on ",
    	                                              sourceSetName,"_",targetSetName,
    	                                              " network",
    	                                              sep=""))
    	    }
    	    targetSetNodeIDs=targetSetNodeIDsList[[systemName]][[targetSetName]]
    	    for(sourceNode in sourceSetNodeIDs) {
    	      for(targetNode in targetSetNodeIDs){
    	        for(weightType in weightTypes){
    	          if(verbose && verboseLevel>1){
    	            print(paste("---loading ",weightType," paths for ",
    	                        sourceNode,"_",targetNode,sep=""))
    	          }
    	          if ( ordered_nodes & as.numeric(sourceNode) > as.numeric(targetNode) ) {
    	            dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
    	                                     sep="_"),
									   paste("window",window,sep="_"),
    	                               weightType,
    	                               paste(targetNode,sourceNode,sep="_"),
    	                               suboptFileExt,sep=".")
    	          } else {
    	            dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
    	                                     sep="_"),
									   paste("window",window,sep="_"),
    	                               weightType,
    	                               paste(sourceNode,targetNode,sep="_"),
    	                               suboptFileExt,sep=".")
    	          }
    	          
    	          dataFilePath=paste(baseDir,dataFileName,sep="/")
    	          if(!file.exists(dataFilePath)) {
    	            print(paste("Error:",dataFilePath,"does not exist. The file will be skipped."))
    	          } else {
    	            tempFrame=read_networkView_suboptFile_to_DataFrame(dataFilePath,
    	                                                               resCount)
    	            if (length(tempFrame)>0){
    	              #subopt filenames are based on base 0 indexing (vmd)
    	              #but data is loaded into data frames with base 1 indexing (R)
    	              #so we need to update the source and target nodes which are set using
    	              #base 0 indexing (needed to match the base 0 indexed subopt filenames)
    	              #Note that X and Y node ids are okay since this was alread taken care of
    	              #in the 'read_networkView_suboptFile_to_DataFrame' function	
					  pathCount=read_networkView_suboptFile_pathCount(dataFilePath,verbose=FALSE) 
    	              tempFrame=tempFrame%>%
    	                mutate(SYSTEM=systemName,WINDOW=window,WEIGHT_TYPE=weightType,
    	                       SOURCE_SET=sourceSetName,
    	                       TARGET_SET=targetSetName,
    	                       SOURCE_NODE=sourceNode+1, 
    	                       TARGET_NODE=targetNode+1,
							   PATHS_COUNT=pathCount)%>%
    	                dplyr::select(SYSTEM,WINDOW,WEIGHT_TYPE,SOURCE_SET,TARGET_SET,
    	                              SOURCE_NODE,TARGET_NODE,PATHS_COUNT,
    	                              X,Y,COUNT)
    	              suboptDataFrame=rbind(suboptDataFrame,tempFrame)
    	            } else {
    	              print(paste("ERROR: Data frame could not be obtained for ",
    	                          dataFilePath,
    	                          "Skipping File."))
    	            }
    	          }
    	        }
    	      }
    	    }
    	  }
    	}
	}
  }
  return(suboptDataFrame)
}

load_windowed_networkView_subopt_nodeCounts <- function(suboptDataDir,suboptFileExt,
                                        systemNames,windows,weightTypes,resCountsList,
                                        sourceSetNamesList,sourceSetNodeIDsList,
                                        targetSetNamesList,targetSetNodeIDsList,
                                        verbose=FALSE,verboseLevel=0,ordered_nodes=FALSE) {
  #reads sets of subopt files created by networkView into a joint data frame
  #
  #all subopt data files should be stored in suboptDataDir and named as:
  #	systemName_sourceSetName_targetSetName.window_windowID.weightType.sourceNode_targetNode.suboptFileExt
  # e.g. wt_FKBP_res375.window_1.CORRweight.203_375.out.out
  #
  #Argument specifications:
  ### systemNames: 			a vector of strings representing the names of each system to load
  ###							subopt network data for
  ###
  ### windows:				a list / array of window numbers or names
  ### weightTypes: 			a vector of strings representing the names of each weighting used to
  ###				  			calculate different subopt path networks
  ###
  ### resCountsList:			list containing the number of residues for each system
  ###							i.e. accessed as: resCountsList[[systemName]]
  ###
  ### sourceSetNamesList:		a 'list()' of vectors of source set names organized by system name
  ###							i.e. accessed as: sourceSetNamesList[[systemName]]
  ###
  ### sourceSetNodeIDsList:	a list of lists 'list(list())' of vectors of source node indices
  ###							organized by system name then source set name.
  ###							i.e. accesssed as: sourceSetNodeIDsList[[systemName]][[sourceSetName]]
  ### targetSetNamesList:		same as sourceSetNamesList, but for target sets
  ### 
  ### targetSetNodeIDsList:	same as sourceSetNodeIDsList, but of target sets
  suboptDataFrame=data.frame(SYSTEM=c(),WINDOW=c(),WEIGHT_TYPE=c(),
                             SOURCE_SET=c(),TARGET_SET=c(),
                             SOURCE_NODE=c(),TARGET_NODE=c(),
							 PATHS_COUNT=c(),
                             X=c(),Y=c(),COUNT=c())
  for(systemName in systemNames){
    if(verbose){print(paste("working on system",systemName))}
    sourceSetNames=sourceSetNamesList[[systemName]]
    targetSetNames=targetSetNamesList[[systemName]]
    resCount=resCountsList[[systemName]]
	for(window in windows) {
		if(verbose){print(paste("-working on window",window))}
    	for(sourceSetName in sourceSetNames){
    	  sourceSetNodeIDs=sourceSetNodeIDsList[[systemName]][[sourceSetName]]
    	  for(targetSetName in targetSetNames){
    	    if(verbose && verboseLevel>0){print(paste("--working on ",
    	                                              sourceSetName,"_",targetSetName,
    	                                              " network",
    	                                              sep=""))
    	    }
    	    targetSetNodeIDs=targetSetNodeIDsList[[systemName]][[targetSetName]]
    	    for(sourceNode in sourceSetNodeIDs) {
    	      for(targetNode in targetSetNodeIDs){
    	        for(weightType in weightTypes){
    	          if(verbose && verboseLevel>1){
    	            print(paste("---loading ",weightType," paths for ",
    	                        sourceNode,"_",targetNode,sep=""))
    	          }
    	          if ( ordered_nodes & as.numeric(sourceNode) > as.numeric(targetNode) ) {
    	            dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
    	                                     sep="_"),
									   paste("window",window,sep="_"),
    	                               weightType,
    	                               paste(targetNode,sourceNode,sep="_"),
    	                               suboptFileExt,sep=".")
    	          } else {
    	            dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
    	                                     sep="_"),
									   paste("window",window,sep="_"),
    	                               weightType,
    	                               paste(sourceNode,targetNode,sep="_"),
    	                               suboptFileExt,sep=".")
    	          }
    	          
    	          dataFilePath=paste(baseDir,dataFileName,sep="/")
    	          if(!file.exists(dataFilePath)) {
    	            print(paste("Error:",dataFilePath,"does not exist. The file will be skipped."))
    	          } else {
    	            tempFrame=read_networkView_suboptFile_nodeCounts_to_DataFrame(dataFilePath,
    	                                                               resCount)
    	            if (length(tempFrame)>0){
    	              #subopt filenames are based on base 0 indexing (vmd)
    	              #but data is loaded into data frames with base 1 indexing (R)
    	              #so we need to update the source and target nodes which are set using
    	              #base 0 indexing (needed to match the base 0 indexed subopt filenames)
    	              #Note that X and Y node ids are okay since this was alread taken care of
    	              #in the 'read_networkView_suboptFile_to_DataFrame' function	
					  pathCount=read_networkView_suboptFile_pathCount(dataFilePath,verbose=FALSE) 
    	              tempFrame=tempFrame%>%
    	                mutate(SYSTEM=systemName,WINDOW=window,WEIGHT_TYPE=weightType,
    	                       SOURCE_SET=sourceSetName,
    	                       TARGET_SET=targetSetName,
    	                       SOURCE_NODE=sourceNode+1, 
    	                       TARGET_NODE=targetNode+1,
							   PATHS_COUNT=pathCount)%>%
    	                dplyr::select(SYSTEM,WINDOW,WEIGHT_TYPE,SOURCE_SET,TARGET_SET,
    	                              SOURCE_NODE,TARGET_NODE,PATHS_COUNT,
    	                              ResNum,COUNT)
    	              suboptDataFrame=rbind(suboptDataFrame,tempFrame)
    	            } else {
    	              print(paste("ERROR: Data frame could not be obtained for ",
    	                          dataFilePath,
    	                          "Skipping File."))
    	            }
    	          }
    	        }
    	      }
    	    }
    	  }
    	}
	}
  }
  return(suboptDataFrame)
}


load_networkView_suboptData <- function(suboptDataDir,suboptFileExt,
                                        systemNames,weightTypes,resCountsList,
                                        sourceSetNamesList,sourceSetNodeIDsList,
                                        targetSetNamesList,targetSetNodeIDsList,
                                        verbose=FALSE,verboseLevel=0,ordered_nodes=FALSE) {
  #reads sets of subopt files created by networkView into a joint data frame
  #
  #all subopt data files should be stored in suboptDataDir and named as:
  #	systemName_sourceSetName_targetSetName.weightType.sourceNode_targetNode.suboptFileExt
  #
  #Argument specifications:
  ### systemNames: 			a vector of strings representing the names of each system to load
  ###							subopt network data for
  ###
  ### weightTypes: 			a vector of strings representing the names of each weighting used to
  ###				  			calculate different subopt path networks
  ###
  ### resCountsList:			list containing the number of residues for each system
  ###							i.e. accessed as: resCountsList[[systemName]]
  ###
  ### sourceSetNamesList:		a 'list()' of vectors of source set names organized by system name
  ###							i.e. accessed as: sourceSetNamesList[[systemName]]
  ###
  ### sourceSetNodeIDsList:	a list of lists 'list(list())' of vectors of source node indices
  ###							organized by system name then source set name.
  ###							i.e. accesssed as: sourceSetNodeIDsList[[systemName]][[sourceSetName]]
  ### targetSetNamesList:		same as sourceSetNamesList, but for target sets
  ### 
  ### targetSetNodeIDsList:	same as sourceSetNodeIDsList, but of target sets
  suboptDataFrame=data.frame(SYSTEM=c(),WEIGHT_TYPE=c(),
                             SOURCE_SET=c(),TARGET_SET=c(),
                             SOURCE_NODE=c(),TARGET_NODE=c(),
                             X=c(),Y=c(),COUNT=c())
  for(systemName in systemNames){
    if(verbose){print(paste("working on system",systemName))}
    sourceSetNames=sourceSetNamesList[[systemName]]
    targetSetNames=targetSetNamesList[[systemName]]
    resCount=resCountsList[[systemName]]
    for(sourceSetName in sourceSetNames){
      sourceSetNodeIDs=sourceSetNodeIDsList[[systemName]][[sourceSetName]]
      for(targetSetName in targetSetNames){
        if(verbose && verboseLevel>0){print(paste("-working on ",
                                                  sourceSetName,"_",targetSetName,
                                                  " network",
                                                  sep=""))
        }
        targetSetNodeIDs=targetSetNodeIDsList[[systemName]][[targetSetName]]
        for(sourceNode in sourceSetNodeIDs) {
          for(targetNode in targetSetNodeIDs){
            for(weightType in weightTypes){
              if(verbose && verboseLevel>1){
                print(paste("--loading ",weightType," paths for ",
                            sourceNode,"_",targetNode,sep=""))
              }
              if ( ordered_nodes & as.numeric(sourceNode) > as.numeric(targetNode) ) {
                dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
                                         sep="_"),
                                   weightType,
                                   paste(targetNode,sourceNode,sep="_"),
                                   suboptFileExt,sep=".")
              } else {
                dataFileName=paste(paste(systemName,sourceSetName,targetSetName,
                                         sep="_"),
                                   weightType,
                                   paste(sourceNode,targetNode,sep="_"),
                                   suboptFileExt,sep=".")
              }
              
              dataFilePath=paste(baseDir,dataFileName,sep="/")
              if(!file.exists(dataFilePath)) {
                print(paste("Error:",dataFilePath,"does not exist. The file will be skipped."))
              } else {
                tempFrame=read_networkView_suboptFile_to_DataFrame(dataFilePath,
                                                                   resCount)
                if (length(tempFrame)>0){
                  #subopt filenames are based on base 0 indexing (vmd)
                  #but data is loaded into data frames with base 1 indexing (R)
                  #so we need to update the source and target nodes which are set using
                  #base 0 indexing (needed to match the base 0 indexed subopt filenames)
                  #Note that X and Y node ids are okay since this was alread taken care of
                  #in the 'read_networkView_suboptFile_to_DataFrame' function
                  tempFrame=tempFrame%>%
                    mutate(SYSTEM=systemName,WEIGHT_TYPE=weightType,
                           SOURCE_SET=sourceSetName,
                           TARGET_SET=targetSetName,
                           SOURCE_NODE=sourceNode+1, 
                           TARGET_NODE=targetNode+1)%>%
                    dplyr::select(SYSTEM,WEIGHT_TYPE,SOURCE_SET,TARGET_SET,
                                  SOURCE_NODE,TARGET_NODE,
                                  X,Y,COUNT)
                  suboptDataFrame=rbind(suboptDataFrame,tempFrame)
                } else {
                  print(paste("ERROR: Data frame could not be obtained for ",
                              dataFilePath,
                              "Skipping File."))
                }
              }
            }
          }
        }
      }
    }
  }
  return(suboptDataFrame)
}
summarise_suboptData<-function(suboptData) {
  return(suboptDataFrame %>%
           group_by(SYSTEM,WEIGHT_TYPE,SOURCE_SET,TARGET_SET,X,Y) %>%
           summarise(COUNT_TOTAL=sum(COUNT)))
}
write_suboptUsageSummaries_to_networkView_Matrices<-function(baseDir,outputExt="subopt.dat",resCountsList,
                                                             suboptData,verbose=FALSE){
  suboptDataSummary<-summarise_suboptData(suboptData)
  for(systemName in suboptDataSummary$SYSTEM%>%factor%>%levels){
    resCount=resCountsList[[systemName]]
    for(weightType in suboptDataSummary$WEIGHT_TYPE%>%factor%>%levels){
      for(sourceSet in suboptDataSummary$SOURCE_SET%>%factor%>%levels){
        for(targetSet in suboptDataSummary$TARGET_SET%>%factor%>%levels) {
          outputNetworkName=paste(systemName,sourceSet,targetSet,sep="_")
          if(verbose){print(paste("Writting data for",outputNetworkName,weightType,"network"))}
          flush.console()
          tempFrame=suboptDataSummary %>% 
            filter(SYSTEM==systemName,SOURCE_SET==sourceSet,TARGET_SET==targetSet,
                   COUNT_TOTAL>0) %>%
            mutate(USAGE=COUNT_TOTAL/max(COUNT_TOTAL),
                   USAGEweight=-log(USAGE),
                   USAGEstrength=-1/(-1+log(USAGE)))
          weightFileName=paste(outputNetworkName,weightType,"USAGEweight",outputExt,sep=".")
          strengthFileName=paste(outputNetworkName,weightType,"USAGEstrength",outputExt,sep=".")
          weightFilePath=paste(baseDir,weightFileName,sep="/")
          strengthFilePath=paste(baseDir,strengthFileName,sep="/")
          write_networkDataTable_to_networkView_matrix(weightFilePath,
                                                       networkData=tempFrame,valColumn="USAGEweight",
                                                       nodeCount=resCount,indexOffset=-1)
          write_networkDataTable_to_networkView_matrix(strengthFilePath,
                                                       networkData=tempFrame,valColumn="USAGEstrength",
                                                       nodeCount=resCount,indexOffset=-1)
        }
      }
    }
  }
  if(verbose){print("DONE")}	
}
###################################################################################################
##########################################IGRAPH section ##########################################
#These utilities allow converting data loaded from carma network analysis files into igraph
#data structures and allow generation of current-flow-betweenness based network analysis
#and relevant plots.
#SETUP info
# data files needed:
# contact_correlation_data (see above utilities to load from carma data files)
# res_data (contains info to map node index to residue index / name and coordinates)
#
#Ex1: Generating plot ready igraph data for system "wt" from node 35 to node 204
#     and label relevant mutation locations at residues 206 and 258, with path data for
#     suboptimal current-flow-betweenness paths up to +25% of optimal path distance 
#plot_graph <- gen_btw_subOptPath_graph_from_Data(res_data,contact_correlation_data,"wt",
#                                   srcnodes=c(35),trgnodes=c(204),mResids = c(206,258),
#                                   dist=.25,byRelativeDist=TRUE,verbose=TRUE)
#
#Ex2: use above plot_graph to plot raw correlation colored / weighted graph
#       par(bg="grey50")
#       plot_full_graph_diverge_hcl(plot_graph,vertSizes = 2,
#                                  edgeSizes = 6*abs(E(plot_graph)$CORR),
#                                  edgeValues= abs(E(plot_graph)$CORR),show_legend=TRUE)
#
#Ex3: use plot_graph to plot graph with edges colored by current-flow-betweenness and 
#     weighted by raw correlation value, and vertex size scaled by current-flow-betweenness
#        par(bg="grey50")
#        plot_full_graph_diverge_hcl(plot_graph,
#                            vertSizes = sqrt(scaleFactor*V(plot_graph)$btw),
#                            edgeSizes = 6*abs(E(plot_graph)$CORR),
#                            edgeValues= sqrt(E(plot_graph)$btw),
#                            show_legend=TRUE)
#
#Ex4: use plot_graph to plot the 25% path dialation current-flow-betweenness sub optimal paths
#        par(bg="grey50")
#        plot_subOpt_flow_matlablike2(plot_graph,
#                             srcnodes = srcnodes,trgnodes = trgnodes,
#                             subOpt_path_data = plot_graph$subOpt_path_data)
#############
#full cutoff
v_btw_from_pInv <- function (Linv,Madj,sIds,tIds,vlist=c()) {
  Linv_dims=dim(Linv);
  Madj_dims=dim(Madj);
  #start Debug
  #print("computing betweeness scores for")
  #print("sIds:")
  #print(sIds)
  #print("tIds:")
  #print(tIds)
  #end Debug
  if(Linv_dims[1]!=Linv_dims[2]) {
    print(paste("ERROR! Pseudo inverse matrix is not square!\nDim1=",
                Linv_dims[1],"; Dim2=",Linv_dims[2],sep=""))
    return(c())
  } else if(Madj_dims[1]!=Madj_dims[2]) {
    print(paste("ERROR! Adjacency matrix is not square!\nDim1=",
                Madj_dims[1],"; Dim2=",Madj[2],sep=""))
    return(c())
  } else if(!(Madj_dims[1]==Linv_dims[1] & Madj_dims[2]==Linv_dims[2])) {
    print(paste("ERROR! Pseudo inverse matrix and Adjacency matrix dimensions do not match!\n Dims1=(",
                Linv_dims[1],",",Madj_dims[1],"); Dim2=(",Linv_dims[2],",",Madj_dims[2],")",sep=""))
    return(c())
  } else if(is.null(vlist) | length(vlist)==0){
    veclist=c(1:Madj_dims[1])
    #print(veclist)
  } else if (max(vlist)>Linv_dims[1] | min(vlist)<1) {
    badvIds=which(vlist>Linv_dims | vlist<1)
    badvs=vlist[badvIds]
    print("ERROR! invalid vertices in vlist!")
    print(paste("--vertexIds:",badvIds,sep=""))
    print(paste("--vertexVals:",badvs,sep=""))
    return(c())
  } else {
    veclist=vlist
  }
  btw=mapply(function(i) {
    sum(c(mapply(function(sId) {
      #start Debug
      #print(paste("-computing btw for sId:",sId))
      #end Debug
      sum(c(mapply(function(tId) {
        gis=Linv[i,sId]
        git=Linv[i,tId]
        gVs=Linv[,sId]
        gVt=Linv[,tId]
        gVabs=abs(gis+gVt-git-gVs)
        1/2*sum(Madj[i,]*gVabs)
      },tIds)))
    },sIds)))
  },veclist)/(length(sIds)*length(tIds))
  return(btw)
}
e_btw_from_pInv <- function(Linv,Madj,sIds,tIds,edgeList) {
  Linv_dims=dim(Linv);
  Madj_dims=dim(Madj);
  if(Linv_dims[1]!=Linv_dims[2]) {
    print(paste("ERROR! Pseudo inverse matrix is not square!\nDim1=",
                Linv_dims[1],"; Dim2=",Linv_dims[2],sep=""))
    return(c())
  } else if(Madj_dims[1]!=Madj_dims[2]) {
    print(paste("ERROR! Adjacency matrix is not square!\nDim1=",
                Madj_dims[1],"; Dim2=",Madj[2],sep=""))
    return(c())
  } else if(!(Madj_dims[1]==Linv_dims[1] & Madj_dims[2]==Linv_dims[2])) {
    print(paste("ERROR! Pseudo inverse matrix and Adjacency matrix dimensions do not match!\n Dims1=(",
                Linv_dims[1],",",Madj_dims[1],"); Dim2=(",Linv_dims[2],",",Madj_dims[2],")",sep=""))
    return(c())
  } 
  Ei=edgeList[,1]
  Ej=edgeList[,2]
  Ebtw=mapply(function(i,j) {
    Madj[i,j] *
      sum(mapply(function(src) {
        sum(mapply(function(trg){
          abs(Linv[i,src]+Linv[j,trg]
              -Linv[i,trg]-Linv[j,src])
        },tIds))
      },sIds))
  },Ei,Ej)/(length(sIds)*length(tIds))
  return(Ebtw)
}
v_btw_kEig <- function(Adj,eVals,eVecs,sIds,tIds,vlist) {
  eCs=mapply(function(eVid){
    1/2*1/eVals[eVid]*sum(
      mapply(function(sId){
        sum(mapply(function(tId){
          abs(eVecs[sId,eVid]-eVecs[tId,eVid])
        },tIds))
      },sIds)
    )
  },c(1:length(eVals)))/(length(sIds)*length(tIds))
  Vbtw=mapply(function(vId){
    sum(mapply(function(eVid){
      eVdelta=abs(eVecs[vId,eVid]-eVecs[,eVid])
      Adj[vId,]*eVdelta*eCs[eVid]
    },c(1:length(eVals))))
  },vlist)
  return(Vbtw)
}
e_btw_kEig <- function(Adj,eVals,eVecs,sIds,tIds,edges) {
  Ei=edges[,1]
  Ej=edges[,2]
}
compute_full_betweenness <- function(g,sourceResNums,targetResNums) {
  #requires the "weight" attribute to be set for edges
  #requires the "ResNum" attribute to be set for vertices
  #the sourceResNums and targetResNums should be subsets of
  #the vertices in g
  
  #sanity check to make sure that only source and target nodes that are
  #in the graph g get used
  srcnodes=V(g)$ResNum[which(is.element(V(g)$ResNum,sourceResNums))]
  trgnodes=V(g)$ResNum[which(is.element(V(g)$ResNum,targetResNums))]
  #start Debug
  #print("computing betweenness scores")
  #print("source residue numbers:")
  #print(sourceResNums)
  #print("corresponding to vertices:")
  #print(srcnodes)
  #print("graph resnums:")
  #print(V(g)$ResNum)
  #end Debug
  Lap=graph.laplacian(g)
  Adj=get.adjacency(g,attr="weight")
  LapInv=as.matrix(Lap) %>% ginv()
  Vbtw=v_btw_from_pInv(Linv = LapInv,Madj = Adj,sIds = srcnodes,tIds=trgnodes)
  ##
  ##Edge Betweeness
  Ebtw=e_btw_from_pInv(LapInv,Adj,srcnodes,trgnodes,
                       get.edgelist(g,names=FALSE))
  return(list(Vbtw=Vbtw,Ebtw=Ebtw))
}
#eigen vector betweenness is currently broken!
#compute_kEig_betweenness <- function(g,sourceResNums,targetResNums,numEigs) {
#  srcnodes=V(g)$ResNum[which(is.element(V(g)$ResNum,sourceResNums))]
#  trgnodes=V(g)$ResNum[which(is.element(V(g)$ResNum,targetResNums))]
#  Lap=graph.laplacian(g)
#  Adj=get.adjacency(g,attr="weight")
#  eigData <- slanczos(Lap,k=0,kl=numEigs+1)
#  eigVals <- eigData$values[1:numEigs]
#  eigVecs <- eigData$vectors[,c(1:numEigs)]
#  Vbtw <- v_btw_kEig(Adj,eigVals,eigVecs,srcnodes,trgnodes,c(1:length(V(g))))
#  
#  return(list(Vbtw=Vbtw))
#}

path_to_eids <- function(g,vpath,verbose=FALSE) {
  if (length(vpath)<2 ) {
    if (verbose) {
      print("WARNING: Call to function 'path_to_eids' with path length<2 detected.")
      print("        Check for disjoint nodes in subgraphs!")
      flush.console()
    }
    return(c())
  }
  path_list=c(1:(2*(length(vpath)-1)))
  path_list[c(1:(length(vpath)-1))*2]=vpath[c(2:(length(vpath)))]
  path_list[c(1:(length(vpath)-1))*2-1]=vpath[c(1:(length(vpath)-1))]
  get.edge.ids(g,path_list)
}
path_through_vert <- function(g,src,trg,vert,eWeights=NULL) {
  E(g)$name=c(1:length(E(g)))
  srcSubNodes=c(1:(trg-1),(trg+1):length(V(g)))
  trgSubNodes=c(1:(src-1),(src+1):length(V(g)))
  src_g=induced.subgraph(g,srcSubNodes)
  src_vpath=V(src_g)$name[
    shortest_paths(src_g,
                   from=which(as.numeric(V(src_g)$name)==src),
                   to=which(as.numeric(V(src_g)$name)==vert),
                   weights=eWeights[E(src_g)$name])$vpath[[1]] %>%
      as.numeric()] %>% as.numeric()
  src_vpath_enames=E(g)$name[path_to_eids(g,src_vpath)]%>%as.numeric()
  trgSubNodes=trgSubNodes[which(
   V(g)$name[trgSubNodes]==V(g)$name[trg] |
     V(g)$name[trgSubNodes]==V(g)$name[vert] |
     !is.element(V(g)$name[trgSubNodes],V(g)$name[src_vpath])
  )] #remove nodes in src_vpath from target pathing subgraph
  trg_g=induced.subgraph(g,trgSubNodes)
  trg_g_edges=E(trg_g)[which(!is.element(
    E(trg_g)$name,src_vpath_enames
  ))] #remove edges of src_vpath from target pathing subgraph
  trg_g=subgraph.edges(trg_g,trg_g_edges)
  trg_vpath=V(trg_g)$name[
    shortest_paths(trg_g,
                   from=which(as.numeric(V(trg_g)$name)==vert),
                   to=which(as.numeric(V(trg_g)$name)==trg),
                   weights=eWeights[E(trg_g)$name])$vpath[[1]] %>%
      as.numeric()] %>% as.numeric()
  src_path_eids=path_to_eids(g,src_vpath)
  trg_path_eids=path_to_eids(g,trg_vpath)
  return(c(src_path_eids,trg_path_eids)%>%unique())
}
distance_ellipse <- function(g,src,trg,eWeights=NULL,dist=0,byRelativeDist=FALSE) {
  E(g)$name=c(1:length(E(g)))
  dmin=as.numeric(distances(g,v=c(src),to=c(trg),weights=eWeights))
  if (byRelativeDist)  {
    dmax=dmin*(1+dist)
  } else {
    dmax=dmin+dist
  }
  srcSubNodes=c(1:(trg-1),(trg+1):length(V(g)))
  trgSubNodes=c(1:(src-1),(src+1):length(V(g)))
  src_g=induced.subgraph(g,srcSubNodes)
  trg_g=induced.subgraph(g,trgSubNodes)
  srcDistances=distances(src_g,v=V(src_g)[which(as.numeric(V(src_g)$name)==src)],
                         to=V(src_g)[which(is.element(as.numeric(V(src_g)$name),srcSubNodes))],
                         weights=eWeights[as.numeric(E(src_g)$name)])%>%as.numeric()
  
  pathCheckNodes=V(src_g)$name[which(srcDistances <= dmax)] %>% as.numeric()
  pathCheckNodes=pathCheckNodes[which(pathCheckNodes != src & pathCheckNodes != trg)]
  dvalList=c(1:length(V(g)))*Inf
  dvalList[trg]=dmin
  dvalList[src]=0
  for(vert in pathCheckNodes) {
    vpath=path_through_vert(g,src,trg,vert,eWeights)
    vpathSub_g=subgraph.edges(g,vpath)
    vpathDist=distances(vpathSub_g,
                        v=V(vpathSub_g)[which(as.numeric(V(vpathSub_g)$name)==src)],
                        to=V(vpathSub_g)[which(as.numeric(V(vpathSub_g)$name)==trg)],
                        weights=eWeights[as.numeric(E(vpathSub_g)$name)])%>%as.numeric()
    dvalList[vert]=vpathDist[1]
  }
  tempframe=data.frame(verts=as.numeric(V(g)$name),dists=dvalList) %>%
    filter(is.finite(dists)&(dists<=dmax))
  return(list(
    vertices=tempframe$verts,
    distances=tempframe$dists
  ))
  #trgDistances=distances(trg_g,v=V(trg_g)[which(as.numeric(V(trg_g)$name)==trg)],
  #                       to=V(trg_g)[which(is.element(as.numeric(V(trg_g)$name),trgSubNodes))],
  #                       weights=eWeights[as.numeric(E(trg_g)$name)])%>%as.numeric()
  #srcDistances=c(srcDistances[c(1:(trg-1))],dmin,srcDistances[c((trg):length(srcDistances))])
  #trgDistances=c(trgDistances[c(1:(src-1))],dmin,trgDistances[c((src):length(trgDistances))])
  #dsum=srcDistances+trgDistances
  #return(list(
  #  vertices=as.numeric(which(dsum<=dmax & is.finite(dsum))),
  #  distances=as.numeric(dsum[which(dsum<=dmax & is.finite(dsum))]))
  #)
}
distance_ellipse_paths_eids <- function(g,src,trg,
                                        eWeights=NULL,dist=0,
                                        byRelativeDist=FALSE) {
  ellipseVerts=distance_ellipse(g,src,trg,eWeights,dist,byRelativeDist)$vertices
  
  E(g)$name=c(1:length(E(g)))
  
  eid_list=c()
  min_vpath=g$name[
    shortest_paths(g,
                   from=src,
                   to=c(trg),
                   weights=eWeights[E(g)$name])$vpath[[1]] %>%
      as.numeric()] %>% as.numeric()
  eid_list=c(eid_list,path_to_eids(g,min_vpath))
  if (length(ellipseVerts > 0)) {
    for (vert in ellipseVerts[which(!is.element(ellipseVerts,c(src,trg)))]) {
      eid_list=c(eid_list,path_through_vert(g,src,trg,vert,eWeights))
    }
    eid_list=c(eid_list,
             shortest_paths(g,from=src,to=trg,weights=eWeights)$vpath[[1]]%>%
               as.numeric()%>%
               path_to_eids(g,.))%>%
      unique()
  }
  return(eid_list)
}
subotimal_paths <- function(g,srcNodes,trgNodes,eWeights=NULL,dist=0,
                            byRelativeDist=FALSE,outputEdges=FALSE,verbose=FALSE) {
  outList=list()
  #these nested for loops are inefficient for large src and trg lists
  count=0
  max_count=length(srcNodes)*length(trgNodes)
  if (max_count > 100) {divCount=100} else if (max_count > 10) {divCount=10} else {divCount=5}
  count_tic=ceiling(max_count/divCount)
  if(verbose){cat("Computing Suboptimal Paths:\n|0%");flush.console()}
  for (src in srcNodes) {
    for (trg in trgNodes) {
      if(verbose){
        if( (count_tic > 1) & (count%%count_tic==(count_tic-1))){
          cat(paste(" ...",
                    sprintf("%3.0f",(1+floor(count/count_tic))*(count_tic/max_count)*100),
                    "%",
                    sep=""))
          if (floor(count/count_tic)%%10==9) {cat("\n   ")}
          flush.console()
        }
        count=count+1
      }
      outList[[paste(src,trg,sep="_")]]=distance_ellipse(g,src,trg,eWeights=eWeights,dist=dist,
                                                         byRelativeDist=byRelativeDist)
      if(outputEdges){
        outList[[paste(src,trg,sep="_")]]$eids=distance_ellipse_paths_eids(
          g,src,trg,eWeights=eWeights,dist=dist,byRelativeDist=byRelativeDist
        )
      }
    }
  }
  if(verbose){cat("|\n");flush.console()}
  return(outList)
}
grid_interp_values <- function(values,gridcount) {
  return((data.frame(X=values) %>% 
            mutate(Y=cut(X,breaks=gridcount,labels=FALSE)))$Y)
}
contact_corr_graph_from_Data <- function(gResDat,corr_data,system,
                                         layoutOrientation="Z") {
  g_contact_corr_graph <- corr_data %>% filter(SYSTEM==system,abs(CORR)>0) %>%
    mutate(weight=1/abs(CORR),
           X_rX=gResDat$rX[X],X_rY=gResDat$rY[X],X_rZ=gResDat$rZ[X],
           Y_rX=gResDat$rX[Y],Y_rY=gResDat$rY[Y],Y_rZ=gResDat$rZ[Y],
           rDIST=sqrt((X_rX - Y_rX)^2 + (X_rY - Y_rY)^2 + (X_rZ-Y_rZ)^2)) %>%
    dplyr::select(X,Y,weight,CORR,CONTACT,rDIST) %>% as.data.frame() %>%
    graph_from_data_frame(directed=FALSE,vertices=gResDat)
  V(g_contact_corr_graph)$ResNum=gResDat$ResNum
  
  if(layoutOrientation=="X") {
    g_contact_corr_graph$layout = 
      data.frame(rX=V(g_contact_corr_graph)$rY,
                 rY=V(g_contact_corr_graph)$rZ)%>%
      data.matrix()
  } else if (layoutOrientation=="Y") {
    g_contact_corr_graph$layout = 
      data.frame(rX=V(g_contact_corr_graph)$rX,
                 rY=V(g_contact_corr_graph)$rZ)%>%
      data.matrix()
  } else {
    g_contact_corr_graph$layout = 
      data.frame(rX=V(g_contact_corr_graph)$rX,
                 rY=V(g_contact_corr_graph)$rY)%>%
      data.matrix()
  }
  
  return(g_contact_corr_graph)
}
corr_flow_btw_graph_from_Data <- function(gResDat,graph_data,
                                         layoutOrientation="Z") {
  g_corr_flow_btw_graph <- graph_data %>%
    mutate(X=from,Y=to) %>%
    mutate(X_rX=gResDat$rX[X],X_rY=gResDat$rY[X],X_rZ=gResDat$rZ[X],
           Y_rX=gResDat$rX[Y],Y_rY=gResDat$rY[Y],Y_rZ=gResDat$rZ[Y],
           rDIST=sqrt((X_rX - Y_rX)^2 + (X_rY - Y_rY)^2 + (X_rZ-Y_rZ)^2)) %>%
    dplyr::select(X,Y,name,weight,CORR,CONTACT,rDIST,btw) %>% as.data.frame() %>%
    graph_from_data_frame(directed=FALSE,vertices=gResDat)
  V(g_corr_flow_btw_graph)$ResNum=gResDat$ResNum
  
  if(layoutOrientation=="X") {
    g_corr_flow_btw_graph$layout = 
      data.frame(rX=V(g_corr_flow_btw_graph)$rY,
                 rY=V(g_corr_flow_btw_graph)$rZ)%>%
      data.matrix()
  } else if (layoutOrientation=="Y") {
    g_corr_flow_btw_graph$layout = 
      data.frame(rX=V(g_corr_flow_btw_graph)$rX,
                 rY=V(g_corr_flow_btw_graph)$rZ)%>%
      data.matrix()
  } else {
    g_corr_flow_btw_graph$layout = 
      data.frame(rX=V(g_corr_flow_btw_graph)$rX,
                 rY=V(g_corr_flow_btw_graph)$rY)%>%
      data.matrix()
  }
  
  return(g_corr_flow_btw_graph)
}
add_src_trg_labeling_to_graph <- function(g,gResDat=as_data_frame(g,what="vertices"),
                                          srcnodes,trgnodes,mResids=c(),
                                          srcColor="orange",
                                          trgColor="red",
                                          mColor="darkviolet") {
  color_list=rep("black",length(gResDat$RESID))
  color_list[which(is.element(gResDat$RESID,mResids))]=mColor
  color_list[which(is.element(gResDat$ResNum,srcnodes))]=srcColor
  color_list[which(is.element(gResDat$ResNum,trgnodes))]=trgColor
  
  label_list=rep("",length(gResDat$RESID))
  #label_list[which(is.element(gResDat$RESID,srcnodes))]=
  #  gResDat$Residue[which(is.element(gResDat$RESID,srcnodes))]
  label_list[which(is.element(gResDat$ResNum,trgnodes))]=
    gResDat$Residue[which(is.element(gResDat$ResNum,trgnodes))]
  label_list[which(is.element(gResDat$RESID,mResids))]=
    gResDat$Residue[which(is.element(gResDat$RESID,mResids))]
  V(g)$color="black"
  V(g)$color=color_list
  V(g)$label=label_list
  V(g)$label.dist=.25
  V(g)$label.font=2
  V(g)$size=3
  V(g)$frame.color=V(g)$color
  V(g)$label.color=V(g)$color
  return(g)
}
gen_plot_graph_from_Data <- function(gResDat,corr_data,system,
                                     srcnodes=c(),trgnodes=c(),mResids=c(),
                                     srcColor="orange",
                                     trgColor="red",
                                     mColor="green") {
  
  g_contact_corr_graph <- contact_corr_graph_from_Data(gResDat,corr_data,system)
  
  add_src_trg_labeling_to_graph(g_contact_corr_graph,gResDat,srcnodes,trgnodes,mResids)
  
  return(g_contact_corr_graph %>% simplify(edge.attr.comb="mean"))
}
gen_plot_graph_from_contact_graph <- function(g,srcnodes=c(),trgnodes=c(),mResids=c(),
                                              srcColor="orange",
                                              trgColor="red",
                                              mColor="darkviolet") {
  gResDat <- as_data_frame(g,what="vertices")
  g <- add_src_trg_labeling_to_graph(g,gResDat,srcnodes,trgnodes,mResids)
  return(g %>% simplify(edge.attr.comb="mean"))
}
plot_full_graph_hcl <- function(g,vertSizes,edgeSizes,edgeValues,
                                show_legend=FALSE,showLabels=TRUE) {
  V(g)$size=vertSizes
  E(g)$width=edgeSizes
  edgeColorInds=grid_interp_values(edgeValues,length(edgeValues))
  edgeColors=sequential_hcl(length(edgeValues))[edgeColorInds]
  E(g)$values=edgeValues
  E(g)$color=edgeColors
  
  g_frames=as_data_frame(g,what="both")
  
  g_frames$edges <- g_frames$edges %>%
    arrange(desc(edgeValues))
  
  g_plot <- graph_from_data_frame(g_frames$edges,directed=FALSE,vertices=g_frames$vertices)
  g_plot$layout=g$layout
  
  if(showLabels==FALSE) {
    V(g_plot)$label=""
  }
  
  plot.igraph(g_plot)
  
  if (show_legend) {
    legend.gradient(cbind(x =c(0.85,0.9,0.9,0.85)+.25, y =c(1.0,1.0,0.8,0.8)*8-7.5),
                    cols=sequential_hcl(length(edgeValues)),
                    limits=round(c(-min(edgeValues),-max(edgeValues)),3),
                    title="")
  }
}
plot_full_graph_heat_hcl <- function(g,vertSizes,edgeSizes,edgeValues,
                                     show_legend=FALSE,showLabels=TRUE) {
  V(g)$size=vertSizes
  E(g)$width=edgeSizes
  edgeColorInds=grid_interp_values(edgeValues,length(edgeValues))
  edgeColors=heat_hcl(length(edgeValues),c=c(80,30),l=c(30,90),power=c(1.5,1.5))[edgeColorInds]
  E(g)$values=edgeValues
  E(g)$color=edgeColors
  
  g_frames=as_data_frame(g,what="both")
  
  g_frames$edges <- g_frames$edges %>%
    arrange(desc(edgeValues))
  
  g_plot <- graph_from_data_frame(g_frames$edges,directed=FALSE,vertices=g_frames$vertices)
  g_plot$layout=g$layout
  
  if(showLabels==FALSE) {
    V(g_plot)$label=""
  }
  
  plot.igraph(g_plot)

  if (show_legend) {
    legend.gradient(cbind(x =c(0.85,0.9,0.9,0.85)+.25, y =c(1.0,1.0,0.8,0.8)*8-7.5),
                    cols=heat_hcl(length(edgeValues),c=c(80,30),l=c(30,90),power=c(1.5,1.5)),
                    limits=round(c(-min(edgeValues),-max(edgeValues)),3),
                    title="")
  }
}
plot_full_graph_diverge_hcl <- function(g,vertSizes,edgeSizes,edgeValues,
                                        show_legend=FALSE,showLabels=TRUE) {
  V(g)$size=vertSizes
  E(g)$width=edgeSizes
  edgeColorInds=grid_interp_values(edgeValues,length(edgeValues))
  edgeColors=diverge_hcl(length(edgeValues),c=c(80,30),l=c(30,90),power=c(1.5,1.5))[edgeColorInds]
  E(g)$values=edgeValues
  E(g)$color=edgeColors
  
  g_frames=as_data_frame(g,what="both")
  
  g_frames$edges <- g_frames$edges %>%
    arrange(edgeValues)
  
  g_plot <- graph_from_data_frame(g_frames$edges,directed=FALSE,vertices=g_frames$vertices)
  g_plot$layout=g$layout
  
  if(showLabels==FALSE) {
    V(g_plot)$label=""
  }
  
  plot.igraph(g_plot)
  
  if (show_legend) {
    legend.gradient(cbind(x =c(0.85,0.9,0.9,0.85)+.25, y =c(1.0,1.0,0.8,0.8)*8-7.5),
                    cols=diverge_hcl(length(edgeValues),c=c(80,30),l=c(30,90),power=c(1.5,1.5)),
                    limits=round(c(min(edgeValues),max(edgeValues)),3),
                    title="")
  }
}
plot_full_graph_matlablike2 <- function(g,vertSizes,edgeSizes,edgeValues,
                                        show_legend=FALSE,showLabels=TRUE) {
  V(g)$size=vertSizes
  E(g)$width=edgeSizes
  edgeColorInds=grid_interp_values(edgeValues,length(edgeValues))
  edgeColors=colorRamps::matlab.like2(length(edgeValues))[edgeColorInds]
  E(g)$values=edgeValues
  E(g)$color=edgeColors
  
  g_frames=as_data_frame(g,what="both")
  
  g_frames$edges <- g_frames$edges %>%
    arrange(edgeValues)
  
  g_plot <- graph_from_data_frame(g_frames$edges,directed=FALSE,vertices=g_frames$vertices)
  g_plot$layout=g$layout
  
  if(showLabels==FALSE) {
    V(g_plot)$label=""
  }
  
  plot.igraph(g_plot)
  
  if (show_legend) {
    legend.gradient(cbind(x =c(0.85,0.9,0.9,0.85)+.25, y =c(1.0,1.0,0.8,0.8)*8-7.5),
                    cols=colorRamps::matlab.like2(length(edgeValues)),
                    limits=round(c(min(edgeValues),max(edgeValues)),3),
                    title="")
  }
}
gen_subOpt_flow_path_data <- function(g,srcnodes,trgnodes,
                                      Vbtw=V(g)$btw,Ebtw=E(g)$btw,
                                      dist=0,byRelativeDist=FALSE,verbose=FALSE) {
  ellipseVerts=c()
  ellipseEids=c()
  for(pair in subotimal_paths(g,srcnodes,trgnodes,
                              eWeights = 1/abs(Ebtw),dist=dist,
                              byRelativeDist = byRelativeDist,outputEdges = TRUE,
                              verbose=verbose)) {
    #ellipseVerts=c(ellipseVerts,pair$vertices)%>%unique()%>%sort()
    if (length(which(is.finite(pair$distances)))>0) {
      ellipseEids=c(ellipseEids,pair$eids[which(is.finite(pair$distances))])%>%unique()
    }
  }
  E(g)$name=c(1:length(E(g)))
  ellipseEids=ellipseEids[which(is.finite(ellipseEids))]
  plot_edge_ids=ellipseEids
  ellipseVerts=c(get.edges(g,ellipseEids))%>%unique()
  ellipseVerts=ellipseVerts[which(is.finite(ellipseVerts))]
  #print(ellipseVerts)
  plot_vert_ids=which(
    is.element(V(g)$name,
               V(induced_subgraph(g,ellipseVerts))$name))
  plot_vert_ids=plot_vert_ids[which(!is.element(plot_vert_ids,c(srcnodes,trgnodes)))]
  
  return(list(vertIds=plot_vert_ids,edgeIds=plot_edge_ids))
}
plot_subOpt_flow_matlablike2 <- function(g,srcnodes,trgnodes,
                                         subOpt_path_data=g$subOpt_path_data,
                                         Vbtw=V(g)$btw,Ebtw=E(g)$btw,
                                         srcColor="orange",trgColor="red",
                                         pvColor="magenta",scaleFactor=64,
                                         labelPathVertResids=FALSE,
                                         vertResidThreshold=.90,verbose=FALSE,
                                         scaleVerts=TRUE,showLabels=TRUE) {
  plot_vert_ids=subOpt_path_data$vertIds
  plot_edge_ids=subOpt_path_data$edgeIds
  V(g)$color[plot_vert_ids]=pvColor
  V(g)$color[srcnodes]=srcColor
  V(g)$color[trgnodes]=trgColor
  V(g)$size=0
  if(scaleVerts){V(g)$size[plot_vert_ids]=sqrt(scaleFactor*Vbtw[plot_vert_ids])}
  V(g)$size[srcnodes]=sqrt(scaleFactor*.25)
  V(g)$size[trgnodes]=sqrt(scaleFactor*.25)
  #V(g)$label[plot_vert_ids]=wtResDat$RESID[plot_vert_ids]
  E(g)$width=0
  E(g)$width[plot_edge_ids]=sqrt(scaleFactor)*abs(E(g)$CORR[plot_edge_ids])
  edge_colorInds=grid_interp_values(Ebtw%>%sqrt(),length(Ebtw))
  E(g)$color[plot_edge_ids]=(colorRamps::matlab.like2(
    length(Ebtw)))[edge_colorInds[plot_edge_ids]]
  E(g)$value=0
  E(g)$value[plot_edge_ids]=edge_colorInds[plot_edge_ids]
  
  if(labelPathVertResids) {
    threshold=(max(V(g)$btw[plot_vert_ids])-
                 min(V(g)$btw[plot_vert_ids]))*vertResidThreshold+
      min(V(g)$btw[plot_vert_ids])
    label_plot_vert_ids=plot_vert_ids[which(V(g)$btw[plot_vert_ids]>=threshold)]
    V(g)$label[label_plot_vert_ids]=V(g)$RESID[label_plot_vert_ids]
  }
  
  g_frames=as_data_frame(g,what="both")
  
  g_frames$edges <- g_frames$edges %>%
    arrange((value))
  
  g_plot <- graph_from_data_frame(g_frames$edges,directed=FALSE,vertices=g_frames$vertices)
  g_plot$layout=g$layout
  
  par(bg="grey50")
  
  if(showLabels==FALSE) {
    V(g_plot)$label=""
  }
  
  plot.igraph(g_plot,
              vertex.label.dist=.25,
              vertex.label.font=1,
              vertex.label.color=V(g_plot)$color %>% as.character(),
              vertex.color=V(g_plot)$color %>% as.character(),
              vertex.frame.color=V(g_plot)$color %>% as.character())
  
  legend.gradient(cbind(x =c(0.85,0.9,0.9,0.85)+.25, y =c(1.0,1.0,0.8,0.8)*8-7.5),
                  cols=colorRamps::matlab.like2(length(Ebtw)),
                  limits=round(c(min(Ebtw),max(Ebtw)),3),
                  title="Normalized\nFlow")
}
gen_btw_subOptPath_graph_from_Data <- function(gResData,corr_data,system,
                                               srcnodes,trgnodes,mResids=c(),
                                               srcColor="orange",trgColor="red",
                                               mColor="darkviolet",
                                               dist=0,byRelativeDist=FALSE,
                                               verbose=FALSE) {
  if(verbose){print("Generating contact correlation graph");flush.console()}
  contact_corr_graph <- contact_corr_graph_from_Data(gResData,corr_data,system)
  if(verbose){print("Generating Current-Flow-Betweenness data");flush.console()}
  btw_dat=compute_full_betweenness(contact_corr_graph,srcnodes,trgnodes)
  V(contact_corr_graph)$btw=btw_dat$Vbtw
  E(contact_corr_graph)$btw=btw_dat$Ebtw
  rm(btw_dat)
  if(verbose) {print("Constructing plot graph from contact graph and betweenness data");flush.console()}
  plot_graph <- gen_plot_graph_from_contact_graph(contact_corr_graph,
                                                  srcnodes,trgnodes,mResids)
  if(verbose) {print("Computing flow path data");flush.console()}
  subOpt_path_data=gen_subOpt_flow_path_data(plot_graph,srcnodes,trgnodes,
                                             dist=dist,byRelativeDist = byRelativeDist,
                                             verbose=verbose)
  plot_graph$subOpt_path_data=subOpt_path_data
  rm(subOpt_path_data)
  if(verbose) {print("Done!");flush.console()}
  return(plot_graph)
}
gen_subOpt_graph_info_Frames <- function(g,subOptData=g$subOpt_path_data)
  {
  vertIds=subOptData$vertIds
  edgeIds=subOptData$edgeIds
  vertInfo=as_data_frame(g,what="vertices")[
    vertIds,c("name","RESID","RESNAME","Residue","rX","rY","rZ","ResNum","btw")]
  edgeInfo=as_data_frame(g,what="edges")[
    edgeIds,c("from","to","weight","CORR","CONTACT","rDIST","btw")]
  return(list(vertInfo=vertInfo,edgeInfo=edgeInfo))
}

gen_btw_graph_from_Data <- function(gResData,corr_data,system,
                                    srcnodes,trgnodes,mResids=c(),
                                    srcColor="orange",trgColor="red",
                                    mColor="darkviolet",
                                    verbose=FALSE) {
  if(verbose){print("Generating contact correlation graph");flush.console()}
  contact_corr_graph <- contact_corr_graph_from_Data(gResData,corr_data,system)
  if(verbose){print("Generating Current-Flow-Betweenness data");flush.console()}
  btw_dat=compute_full_betweenness(contact_corr_graph,srcnodes,trgnodes)
  V(contact_corr_graph)$btw=btw_dat$Vbtw
  E(contact_corr_graph)$btw=btw_dat$Ebtw
  rm(btw_dat)
  if(verbose) {print("Constructing plot graph from contact graph and betweenness data");flush.console()}
  plot_graph <- gen_plot_graph_from_contact_graph(contact_corr_graph,
                                                  srcnodes,trgnodes,mResids)
  if(verbose) {print("Done!");flush.console()}
  return(plot_graph)
}

do_windowed_betweenness_analysis <- function(windowed_corr_data,resData,
                                             sourceSetsList,targetSetsList,
                                             sourceNodeList,targetNodeList,
                                             nRes=max(c(windowed_corr_data$X,windowed_corr_data$Y)),
                                             windows=windowed_corr_data$WINDOW %>%
                                               unique %>% as.character,
                                             systems=windowed_corr_data$SYSTEM %>%
                                               unique %> %as.character,
                                             includeGraphList=FALSE,mResids=c()) {
  if (includeGraphList==TRUE){
    graphList=list()
  }
  count=0
  for(system in systems){
    resDat=resData %>% filter(SYSTEM==system) %>%
      dplyr::select(ResNum,RESID,RESNAME,Residue,rX,rY,rZ,SYSTEM)
    for(window in windows){
      contactDat=windowed_corr_data %>%
        filter(SYSTEM==system,WINDOW==window)
      print("------------------------------")
      sources=sourceSetsList[[system]]
      targets=targetSetsList[[system]]
      print(paste("Working on",system,"system"))
      print("-source sets:")
      print(sources)
      print("-target sets:")
      for(source in sources) {
        srcnodes=sourceNodeList[[source]]
        for(target in targets) {
          trgnodes=targetNodeList[[target]]
          print("--Source nodes:")
          print(srcnodes)
          print("--Target nodes:")
          print(trgnodes)
          flush.console()
          temp_graph=gen_btw_graph_from_Data(resDat,contactDat,system,
                                             srcnodes = srcnodes,trgnodes = trgnodes,
                                             mResids = mResids,verbose=TRUE)
          if(includeGraphList==TRUE){
            graphList[paste(paste(system,source,target,sep="_"),
                            paste("window",window,sep="_"),sep=".")]=temp_graph
          }
          temp_frame <- as_data_frame(temp_graph) %>% 
            rename(X=from,Y=to)%>%
            mutate(X_Y=paste(X,Y,sep="_")) %>%
            mutate(X=X%>%as.numeric,Y=Y%>%as.numeric,
                   SYSTEM=system,WINDOW=window,NETWORK=paste(source,target,sep="_")) %>%
            select(SYSTEM,WINDOW,NETWORK,X,Y,X_Y,CORR,btw)
          if(count==0){
            network_analysis_frame=temp_frame
          } else {
            network_analysis_frame=rbind(network_analysis_frame,temp_frame)
          }
          count=count+1
        }
      }
    }
  }
  if(includeGraphList) {
    return(list(graphList=graphList,dataFrame=network_analysis_frame))
  } else {
    return(network_analysis_frame)
  }
}
write_tcl_load_scripts <- function(systems,sourceSetsList,targetSetsList,
                                   sourceNodeList,targetNodeList,weightTypes,
                                   outdir="",weightdir="",
                                   verbose=FALSE,verboseLevel=0,
                                   loaderSource="load_suboptimal_paths.tcl") {
  for(system in systems) {
    sourceSet = sourceSetsList[[system]]
    targetSet = targetSetsList[[system]]
    for(source in sourceSet) {
      sourceNodes=sourceNodeList[[source]]
      for(target in targetSet) {
        targetNodes=targetNodeList[[target]]
        for(weightType in weightTypes) {
          outFileName=paste(paste(system,source,target,sep="_"),
                            weightType,"subopt.loader.tcl",sep=".")
          outFilePath=paste(outdir,outFileName,sep="/")
          tcl_lines=c(paste("source",loaderSource),
                      paste("set sourceNodeIDs {",
                            paste(sourceNodes,collapse=" "),"}"),
                      paste("set targetNodeIDs {",
                            paste(targetNodes,collaps=" "),"}"),
                      paste("load_suboptimal_paths",outdir,outdir,system,
                            source,target,"$sourceNodeIDs","$targetNodeIDs",
                            weightType))
          writeLines(tcl_lines,outFilePath)
          if(verbose & verboseLevel>1) {
            print(paste(paste(rep(" - ",5),collapse=""),
                        outFileName,
                        paste(rep(" - ",5),collapse="")))
            print(paste("filepath:",outFilePath))
            for(line in tcl_lines){
              print(line)
            }
            print(paste(rep(" - ",10+round(nchar(outFileName)/3,0)+1),collapse=""))
          }
        }
      }
    }
  }
}
run_external_subopt <- function(systems,sourceSetsList,targetSetsList,
                                sourceNodeList,targetNodeList,
                                weightTypes,dilationList,
                                outdir="",weightdir="",suboptexe="subopt",
                                verbose=FALSE,verboseLevel=0){
  for(system in systems) {
    if(verbose){print(paste("Working on",system,"system"))}
    sourceSet=sourceSetsList[[system]]
    targetSet=targetSetsList[[system]]
    for(source in sourceSet) {
      for(target in targetSet) {
        networkName=paste(source,target,sep="_")
        if(verbose & verboseLevel>0){print(paste("-working on",networkName,"network"))}
        for(weightType in weightTypes){
          weightFileName=paste(paste(system,networkName,sep="_"),
                               weightType,"dat",sep=".")
          weightFilePath=paste(weightdir,weightFileName,sep="/")
          if(!file.exists(weightFilePath)){
            print(paste("ERROR: the weightFile",
                        weightFilePath,"does not exist!"))
            next
          }
          sourceNodes=sourceNodeList[[source]]
          targetNodes=targetNodeList[[target]]
          dilation=dilationList[[weightType]]
          for(sourceNode in sourceNodes) {
            for(targetNode in targetNodes) {
              if (as.numeric(sourceNode) <= as.numeric(targetNode)) {
                nodePairName=paste(sourceNode,targetNode,sep="_")
                sourceID=sourceNode
                targetID=targetNode
              } else {
                nodePairName=paste(targetNode,sourceNode,sep="_")
                sourceID=targetNode
                targetID=sourceNode
              }
              if(verbose & verboseLevel>1){print(paste("---running subopt program for the",
                                                       nodePairName,"source-target pair"))}
              outFileName=paste(paste(system,networkName,sep="_"),
                                weightType,nodePairName,
                                "subopt.out",sep=".")
              outFilePath=paste(outdir,outFileName,sep="/")
              if(verbose & verboseLevel>2) {
                print(paste("----weightFilePath:",weightFilePath))
                print(paste("----outFilePath:",outFilePath))
              }
              system2(suboptexe,
                      args=c(weightFilePath,outFilePath,
                             dilation,sourceID,targetID))
            }
          }
        }
      }
    }
  }
}

run_windowed_external_subopt <- function(systems,windows,sourceSetsList,targetSetsList,
                                sourceNodeList,targetNodeList,
                                weightTypes,dilationList,
                                outdir="",weightdir="",suboptexe="subopt",
                                verbose=FALSE,verboseLevel=0){
  for(system in systems) {
    if(verbose){print(paste("Working on",system,"system"))}
    sourceSet=sourceSetsList[[system]]
    targetSet=targetSetsList[[system]]
    for(source in sourceSet) {
      for(target in targetSet) {
        networkName=paste(source,target,sep="_")
		if(verbose & verboseLevel>0){print(paste("-working on",networkName,"network"))}
		for(window in windows) {
        	if(verbose & verboseLevel>0){print(paste("--working on window",window))}
        	for(weightType in weightTypes){
        	  weightFileName=paste(paste(system,networkName,sep="_"),
								   paste("window",window,sep="_"),
        	                       weightType,"dat",sep=".")
        	  weightFilePath=paste(weightdir,weightFileName,sep="/")
        	  if(!file.exists(weightFilePath)){
        	    print(paste("ERROR: the weightFile",
        	                weightFilePath,"does not exist!"))
        	    next
        	  }
        	  sourceNodes=sourceNodeList[[source]]
        	  targetNodes=targetNodeList[[target]]
        	  dilation=dilationList[[weightType]]
        	  for(sourceNode in sourceNodes) {
        	    for(targetNode in targetNodes) {
        	      if (as.numeric(sourceNode) <= as.numeric(targetNode)) {
        	        nodePairName=paste(sourceNode,targetNode,sep="_")
        	        sourceID=sourceNode
        	        targetID=targetNode
        	      } else {
        	        nodePairName=paste(targetNode,sourceNode,sep="_")
        	        sourceID=targetNode
        	        targetID=sourceNode
        	      }
        	      if(verbose & verboseLevel>1){print(paste("---running subopt program for the",
        	                                               nodePairName,"source-target pair"))}
        	      outFileName=paste(paste(system,networkName,sep="_"),
									paste("window",window,sep="_"),
        	                        weightType,nodePairName,
        	                        "subopt.out",sep=".")
        	      outFilePath=paste(outdir,outFileName,sep="/")
        	      if(verbose & verboseLevel>2) {
        	        print(paste("----weightFilePath:",weightFilePath))
        	        print(paste("----outFilePath:",outFilePath))
        	      }
        	      system2(suboptexe,
        	              args=c(weightFilePath,outFilePath,
        	                     dilation,sourceID,targetID))
        	    }
        	  }
        	}
		}
      }
    }
  }
}


do_windowed_network_aov <- function(network_analysis_data,
                                    outfilebase="all_aov_summary",
                                    chunkDir="all_aov_chunks",
                                    chunksize=30,verbose=FALSE,verboseLevel=0) {
  #prep data from network_analysis_frame for aov
  #need to make sure every edge present has an entry for every window
  #missing entries are set to 0 (since a missing entry most likely means lack of contact)
  if(verbose){
    print("Preparing data for analysis")
  }
  aov_data <- network_analysis_data %>%
    melt(id.vars = c("SYSTEM","WINDOW","NETWORK","X","Y","X_Y"),
         variable.name = "MEASUREMENT", value.name="VALUE") %>%
    unite(SYSTEM.NETWORK.MEASUREMENT.WINDOW,SYSTEM,NETWORK,MEASUREMENT,WINDOW,sep=".") %>%
    tidyr::spread(key=SYSTEM.NETWORK.MEASUREMENT.WINDOW,value=VALUE,fill=0) %>%
    melt(id.vars=c("X","Y","X_Y"),
         variable.name="SYSTEM.NETWORK.MEASUREMENT.WINDOW",value.name="VALUE") %>%
    as.data.frame %>%
    tidyr::separate(col=SYSTEM.NETWORK.MEASUREMENT.WINDOW,
                    into=c("SYSTEM","NETWORK","MEASUREMENT","WINDOW"),sep="[.]") %>%
    select(SYSTEM,WINDOW,NETWORK,X,Y,X_Y,MEASUREMENT,VALUE)
  
  if(verbose){
    if(verboseLevel>0){
      print(paste("-Measurement types to be analyzed:",
                  aov_data$MEASUREMENT%>% unique %>% as.character))
    }
    print("Computing windowed averages and pValues")
  }
  #aov_data %>% head
  count=0
  for(measurement in aov_data$MEASUREMENT %>% unique %>% as.character) {
    if(verbose & verboseLevel > 0) {
      paste("working on",measurement) %>% print
    }
    temp_dat <- aov_data %>% filter(MEASUREMENT==measurement) %>%
      group_by(SYSTEM,NETWORK,X,Y,X_Y,MEASUREMENT)
    
    #temp_dat %>% head %>% print
    
    avg_dat <- temp_dat %>%
      summarise(VALUE=mean(VALUE)) %>%
      as.data.frame()
    
    #avg_dat %>% head %>% print
    #temp_dat %>% head %>% print
    
    p_dat <- temp_dat %>% 
      summarise(VALUE=t.test(VALUE)$p.value) %>% as.data.frame %>% 
      #        head %>% print
      mutate(MEASUREMENT=paste(MEASUREMENT,"pValue",sep="."))
    
    if(verbose & verboseLevel > 1){
      print(paste("-",measurement,"pvalue frame preview:"))
      p_dat %>% head %>% print
    }
    if (count==0) {
      aov_summary <- rbind(avg_dat,p_dat) %>% as.data.frame()
    } else {
      aov_summary <- rbind(rbind(aov_summary,avg_dat),p_dat)
    }
    #aov_summary %>% head %>% print
    #typeof(aov_summary) %>% print
    count=count+1
  }
  
  
  
  contacts <- (aov_summary %>%
                 filter(grepl("pValue",MEASUREMENT),is.finite(VALUE)))$X_Y %>% 
    unique %>% as.character
  
  measurements <- (aov_data$MEASUREMENT%>%unique%>%as.character)
  
  if(verbose){
    print("Running anova for each system and contact pair")
    if(verboseLevel>0){
      cat(paste( (measurements%>%length)*(contacts%>%length) / chunksize,"chunks"))
      }
  }
  syscount=aov_data$SYSTEM%>%unique%>%length
  if(verbose & verboseLevel>0){
    cat("\nsyscount = ")
    cat(syscount)
    cat("\n---run summary---\n")
    cat("\n. = good (all systems present)\n")
    cat("o = incomplete (at least 2, but not all, systems present)\n")
    cat("X = uncomputable (fewer than 2 systems preset)\n")
    cat("\n working on chunk \n")
  }
  count=0
  cnum=1
  if(verbose & verboseLevel>0){
    cat(sprintf("%5s ",cnum))
  }
  networks = aov_data$NETWORK %>% unique %>% as.character()
  tempchunk=data.frame(X_Y=c(),X=c(),Y=c(),
                       NETWORK=c(),MEASUREMENT_SYSTEM=c(),MEASURMENT_TYPE=c(),VALUE=c())
  
  outfile=paste(chunkDir,"/",outfilebase,".chunk_",cnum,".dat",sep="")
  
  aov_summary <- aov_summary %>% 
    rename(MEASUREMENT_TYPE=MEASUREMENT,MEASUREMENT_SYSTEM=SYSTEM) %>%
    select(X_Y,X,Y,NETWORK,MEASUREMENT_SYSTEM,MEASUREMENT_TYPE,VALUE)
  aov_data <- aov_data %>% 
    rename(MEASUREMENT_TYPE=MEASUREMENT,MEASUREMENT_SYSTEM=SYSTEM) %>%
    select(X_Y,X,Y,WINDOW,NETWORK,MEASUREMENT_SYSTEM,MEASUREMENT_TYPE,VALUE)
  for(network in networks) {
    if(verbose & verboseLevel==0) {print(paste("-working on",network,"network"))}
    for(measurement in measurements) {
      if(verbose & verboseLevel==0){print(paste("--working on",measurment,"measurements"))}
      for(contact in contacts) {
        if(verbose & verboseLevel==0){print(paste("--working on contact pair",contact))}
        test_dat <- aov_data %>% 
          filter(NETWORK==network,MEASUREMENT_TYPE==measurement,X_Y==contact,
                 is.finite(VALUE))
        sys_count <- test_dat$MEASUREMENT_SYSTEM %>% unique %>% length
        temp_split=strsplit(x=contact,split="_")[[1]]
        temp_x=temp_split[1]
        temp_y=temp_split[2]
        if((sys_count) > 1) {
          if(verbose & verboseLevel>1)  {
            if (sys_count == syscount) {
              cat(sprintf("%1s ","."))
            } else {
              cat(sprintf("%1s ","o"))
            }
          }
          test_aov <- test_dat %>% aov(VALUE~MEASUREMENT_SYSTEM,data=.)
          test_tuk <- (test_aov %>% TukeyHSD)$MEASUREMENT_SYSTEM %>% 
            as.data.frame %>%
            tibble::rownames_to_column(var="MEASUREMENT_SYSTEM")
          temp_dat <- test_tuk %>%
            mutate(X_Y=contact,X=temp_x,Y=temp_y,NETWORK=network,
                   MEASUREMENT_STAT=measurement) %>%
            melt(id.vars=c("X_Y","X","Y","NETWORK",
                           "MEASUREMENT_SYSTEM","MEASUREMENT_STAT"),
                 variable.name="MEASUREMENT",
                 value.name="VALUE")%>%
            mutate(MEASUREMENT=gsub(" ","_",MEASUREMENT)) %>%
            mutate(MEASUREMENT_TYPE=paste(MEASUREMENT_STAT,MEASUREMENT,sep=".")) %>%
            select(X_Y,X,Y,NETWORK,MEASUREMENT_SYSTEM,MEASUREMENT_TYPE,VALUE)
          tempchunk <- rbind(tempchunk,temp_dat)
          aov_summary <- rbind(aov_summary,temp_dat)
        } else {
          if(verbose & verboseLevel > 1) {cat(sprintf("%1s ","X"))}
        }
        count=count+1
        if(count > chunksize) {
          if(verbose & verboseLevel > 1) {
            cat("\n")
          }
          count=0
          cnum=cnum+1
          if(verbose & verboseLevel > 0) {
            cat(sprintf("%5s ",cnum))
          }
          write.table(tempchunk,outfile,quote=FALSE,row.names=FALSE)
          tempchunk=data.frame(X_Y=c(),X=c(),Y=c(),
                               MEASUREMENT_SYSTEM=c(),MEASURMENT_TYPE=c(),VALUE=c())
          outfile=paste(chunkDir,"/",outfilebase,".chunk_",cnum,".dat",sep="")
        }
      }
    }
  }
  
  return(aov_summary)
}
###################################################################################################
