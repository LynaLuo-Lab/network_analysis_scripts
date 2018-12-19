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

utilities_path="correlation_network_analysis_functions.r"
source(utilities_path)

#variables to control file paths for loading / saving of data
basepath="/home/wesley/alk_network_analysis/windowed_analysis"
baseinfodir="/home/wesley/alk_network_analysis/windowed_analysis"
datafilename="network_varcov_table.dat"
contactmapname="contact.contact.table.dat"
coordfilename="carma.average.pdb"
resinfofilename="resinfo.cpptraj.dat"
atominfofilename="atominfo.cpptraj.dat"
system_list=c("wt","R206H","Q207E","Q207D")
window_list=c(0:6) #only needed when generating graphs from windowed carma data
#resid_map=c(172:499,2:108)

#Variables to concrol generating igraph network graphs and plots
scaleFactor=64
systems=system_list #in case you want to only plot some of the systems and not others
sources=c("FKBP") #c("res206","res207","FKBP")
targets=c("res375")
sourcenode_list=list()
sourcenode_list[["FKBP"]]=c(353,363,364,365,368,369,371,373,374,375,376,379,380,381,
                            382,383,384,386,408,409,411,412,413,414,415,416,417,418,420,424,426,434)
#sourcenode_list[["res206"]]=c(35)
sourcenode_list[["res207"]]=c(36)
targetnode_list=list()
targetnode_list[["res375"]]=c(204)
#mResids=c(206,258,356,328,196,207,197,202)
mResids=c(206,207,375) #c(86,206,207,244,246,354,375) 
dist_list=list()
dist_list[["FKBP_res375"]]=.1
#dist_list[["res206_res375"]]=.25
dist_list[["res207_res375"]]=.25
#ensure equal graph layouts
btw_subOptData_list=list()
btw_graph_list=list()

generate_graphs_from_carma=FALSE
if(generate_graphs_from_carma==TRUE) {
#generate contact correlation data frame
data_list <- load_windowed_correlation_data(basepath=basepath,baseinfodir=baseinfodir,
                                            datafilename=datafilename,
                                            contactmapname=contactmapname,
                                            resinfofilename=resinfofilename,
                                            atominfofilename=atominfofilename,
                                            coordfilename=coordfilename,
                                            system_list=system_list,window_list=window_list,
                                            datafiles_as_arrays=FALSE,
                                            contactfiles_as_arrays=FALSE,
                                            verbose=TRUE)
alk_windowed_corr_data <- data_list[["windowed_corr_data"]]
alk_windowed_contact_map <- data_list[["windowed_contact_data"]]
resData <- data_list[["res_data"]]
#rm(data_list)
alk_contact_corr <- gen_contact_corr_summary(alk_windowed_corr_data,alk_windowed_contact_map,
                                             verbose=TRUE)

contactDat=alk_contact_corr
#generate graphs, compute sub optimal paths from flow-betweenness, and plot the graphs
baseResData=resData %>% filter(SYSTEM=="wt") %>%
  dplyr::select(ResNum,RESID,RESNAME,Residue,rX,rY,rZ,SYSTEM)
for(system in systems) {
  print("------------------------------")
  for(source in sources) {
    srcnodes=sourcenode_list[[source]]
    for(target in targets) {
      trgnodes=targetnode_list[[target]]
      dist=dist_list[[paste(source,target,sep="_")]]
      resDat=resData %>% filter(SYSTEM==system) %>%
        dplyr::select(ResNum,RESID,RESNAME,Residue,rX,rY,rZ,SYSTEM)
      resDat$rX=baseResData$rX
      resDat$rY=baseResData$rY
      resDat$rZ=baseResData$rZ
      #start Debug
      #print(paste("working on source set for",source,"and target set",target))
      #print("source set contents:")
      #print(srcnodes)
      #print("target set contents:")
      #print(trgnodes)
      #end Debug
      print(paste("Generating graph for the",system,source,"to",target,"network",sep=" "))
      plot_graph <- gen_btw_subOptPath_graph_from_Data(resDat,contactDat,system,
                                                       srcnodes = srcnodes,trgnodes = trgnodes,
                                                       mResids = mResids,
                                                       dist=dist,byRelativeDist = TRUE,verbose=TRUE)
      V(plot_graph)$label.cex=3
      print(paste("saving",sprintf("%4.1f%s",dist*100,"%"),
                  "dilation sub-optimal-flow path network",
                  sep=" "))
      btw_subOptData_list[[paste(system,
                                 source,
                                 target,
                                 sep="_")]]=gen_subOpt_graph_info_Frames(plot_graph)
      print(paste("saving igraph data"))
      btw_graph_list[[paste(system,source,target,sep="_")]]=plot_graph
      print("Generating Plots")
      par(bg="grey50")
      plot_full_graph_diverge_hcl(plot_graph,
                                  vertSizes = 2,
                                  edgeSizes = 6*abs(E(plot_graph)$CORR),
                                  edgeValues= abs(E(plot_graph)$CORR),
                                  show_legend=TRUE)
      par(bg="grey50")
      plot_full_graph_diverge_hcl(plot_graph,
                                  vertSizes = sqrt(scaleFactor*V(plot_graph)$btw),
                                  edgeSizes = 6*abs(E(plot_graph)$CORR),
                                  edgeValues= sqrt(E(plot_graph)$btw),
                                  show_legend=TRUE)
      par(bg="grey50")
      plot_subOpt_flow_matlablike2(plot_graph,
                                   srcnodes = srcnodes,trgnodes = trgnodes,
                                   subOpt_path_data = plot_graph$subOpt_path_data)
      print("------------------------------")
    }
  }
}
rm(contactDat)
rm(resDat)

#save the network and sub optimal path subnetwork data to files
##Adjust networkViewHeader to match what was used in the network.config file to generate
##The original network. This should allow the matrices to be imported by the VMD
##networkView plugin.
networkViewHeader="((protein) and (not hydrogen)) and (name CA) )"
for(system in systems) {
  for(source in sources) {
    for(target in targets) {
      networkName=paste(system,source,target,sep="_")
      networkGraph=btw_graph_list[[networkName]]
      subOptData=btw_subOptData_list[[networkName]]
      subOptEdgeData=as.data.frame(subOptData$edgeInfo,stringsAsFactors=TRUE)
      nVerts=length(V(networkGraph))
      degLap=laplacian_matrix(networkGraph,weights=NA)
      subOptEdgeCorrMat=sparseMatrix(i=as.numeric(subOptEdgeData$from),
                                    j=as.numeric(subOptEdgeData$to),
                                    x=as.numeric(subOptEdgeData$CORR),
                                    symmetric = TRUE,dims = c(nVerts,nVerts))
      subOptEdgeData$avgDeg <- mapply(function(x,y){
          sqrt(degLap[x,x]*degLap[y,y])
        },as.numeric(subOptEdgeData$from),as.numeric(subOptEdgeData$to))
      subOptEdgeBtwMat=sparseMatrix(i=as.numeric(subOptEdgeData$from),
                                 j=as.numeric(subOptEdgeData$to),
                                 x=as.numeric(1/(-log(subOptEdgeData$btw)) * subOptEdgeData$avgDeg),
                                 symmetric = TRUE,dims = c(nVerts,nVerts))
      outSubOptBtwPath=paste(basepath,system,
                             paste(networkName,".SubOptBtw.dat",sep=""),sep="/")
      write(networkViewHeader,file=outSubOptBtwPath)
      outSubOptCorrPath=paste(basepath,system,
                              paste(networkName,".SubOptCorr.dat",sep=""),sep="/")
      write(networkViewHeader,file=outSubOptCorrPath)
      print(paste("writting ",outSubOptCorrPath))
      outFile=(outSubOptCorrPath,"a");write.matrix(subOptEdgeCorrMat,outSubOptCorrPath);close(outFile)
      print(paste("writting ",outSubOptBtwPath))
      outFile=(outSubOptBtwPath,"a");write.matrix(subOptEdgeBtwMat,outSubOptBtwPath);close(outFile)
      outSubOptTablePath=paste(basepath,system,
                               paste(networkName,".SubOptTable.dat",sep=""),sep="/")
      print(paste("writting ",outSubOptTablePath))
      subOptEdgeData$name=subOptEdgeData%>%rownames()
      tabCols=subOptEdgeData%>%colnames()
      subOptEdgeData<-subOptEdgeData%>%
        dplyr::select(length(tabCols),c(1:(length(tabCols)-1)))
      write.table(subOptEdgeData,outSubOptTablePath,quote=FALSE,row.names=FALSE)
      networkData=as_data_frame(networkGraph,what="edges")
      networkData$names=networkData%>%rownames()
      tabCols=networkData%>%colnames()
      networkData<-networkData%>%
        dplyr::select(length(tabCols),c(1:(length(tabCols)-1)))
      outNetworkTablePath=paste(basepath,system,
                                paste(networkName,".NetworkTable.dat",sep=""),sep="/")
      print(paste("writting ",outNetworkTablePath))
      write.table(networkData,outNetworkTablePath,quote=FALSE,row.names=FALSE)
    }
  }
}

} else {
#load network and suboptimal path data from files and generate plots
networkDataList<-load_corr_flow_btw_tables(basepath,baseinfodir=basepath,
                          networkDataExt="NetworkTable.dat",
                          subOptDataExt="SubOptTable.dat",
                          resinfofilename="resinfo.cpptraj.dat",
                          atominfofilename="atominfo.cpptraj.dat",
                          coordfilename="carma.average.pdb",
                          system_list=systems,src_list=sources,trg_list=targets,
                          loadSubOptData=TRUE,
                          verbose=TRUE)

print("Plotting graphs of loaded network data")
baseResData=networkDataList[[systems[1]]][["vertexData"]]
for(sys in system_list) {
  sysData=networkDataList[[sys]]
  resDat=sysData[["vertexData"]]
  ###
  #Ensure consistent node layout between systems.
  #Note:this only works if the systems have the same number of nodes in the same order
  resDat$rX=baseResData$rX
  resDat$rY=baseResData$rY
  resDat$rZ=baseResData$rZ
  resDat$name=c(1:length(resDat$rX))
  ###
  print("------------------------------")
  
  for(src in sources){
    for(trg in targets){
      srcnodes=sourcenode_list[[src]]
      trgnodes=targetnode_list[[trg]]
      flowName=paste(src,trg,sep="_")
      flowData=sysData[[flowName]]
      networkData=flowData[["network"]]
      subOptEdges=flowData[["subOptEdges"]]
      print(paste("-Generating graph for the",sys,src,"to",trg,"network",sep=" "))
      dataCols=networkData%>%colnames
      plot_graph=networkData %>%
        dplyr::select(c(5,6,7,14:length(dataCols))) %>%
        graph_from_data_frame(directed=FALSE,
                              vertices=data.frame(
                                name=resDat$name,
                                rX=resDat$rX,
                                rY=resDat$rY#,
                                #rZ=resDat$rZ,
                                #ResNum=resDat$ResNum,
                                #RESNAME=resDat$RESNAME,
                                #RESID=resDat$RESID,
                                #Residue=resDat$Residue
                                ))
      print("--making layout")
      plot_graph$layout = 
        data.frame(rX=V(plot_graph)$rX,
                   rY=V(plot_graph)$rY)%>%
        data.matrix()
      plot_graph <- plot_graph %>%
        add_src_trg_labeling_to_graph(gResDat=resDat,
                                      srcnodes=srcnodes,trgnodes=trgnodes,
                                      mResids=mResids,
                                      srcColor="orange",trgColor="red",
                                      mColor="darkviolet")
      print(paste("-Generating sub optimal flow path subgraph"))
      subOpt_graph=plot_graph %>% subgraph.edges(subOptEdges)
      print("--Generating Plots")
      #par(bg="grey50")
      #print("--Correlation Plot")
      #plot_full_graph_diverge_hcl(plot_graph,
      #                            vertSizes = 2,
      #                            edgeSizes = 6*abs(E(plot_graph)$CORR),
      #                            edgeValues= abs(E(plot_graph)$CORR),
      #                            show_legend=TRUE)
      #par(bg="grey50")
      #print("--Flow-Betweenness Plot")
      #plot_full_graph_diverge_hcl(plot_graph,
      #                            vertSizes = 1,#sqrt(scaleFactor*V(plot_graph)$btw),
      #                            edgeSizes = 6*abs(E(plot_graph)$CORR),
      #                            edgeValues= sqrt(E(plot_graph)$btw),
      #                            show_legend=TRUE)
      par(bg="grey50")
      print("--SubOptimal Path Plot")
      plot_subOpt_flow_matlablike2(plot_graph,
                                   srcnodes = srcnodes,trgnodes = trgnodes,
                                   subOpt_path_data = list(
                                     vertIds=V(subOpt_graph)$name%>%as.numeric(),
                                     edgeIds=E(subOpt_graph)$name%>%as.numeric()
                                     ),scaleVerts=FALSE)
      print("------------------------------")
    }
  }
}
}