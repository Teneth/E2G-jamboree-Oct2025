#### Booleans, pct overlap, min dist, and dist matrix for ChromHMM Cell states - JUSTIN LANGERMAN 2025

### Annotation from 
###  Universal annotation of the human genome through integration of over a thousand epigenomic datasets
# Ha Vu & Jason Ernst 


#### Require data loaded in ###########################################################
# E2G_file_location <- "E2G.file.txt"
# Sample_df <- fread(E2G_file_location)


## Dependencies

library(dplyr)
library(data.table)
library(tidygenomics)  ###Throws warnings on genome_intersect, works as of 2025
library(reshape2)  ## dcast/melt throw depreciated warnings, works as of 2025




# ### Choose assays_calculating wanted. This function can do 5 things:
# 1 "boolean_overlap" ## Calculates the most simply presence or absence of region to state overlap, given as 1 or 0
# 2 "percent_overlap" ## Calculates the % overlap for each region and state, 0 to 1 (100%), also defines majority state
# 2b "majority_state" Calculate winning max % feature
# 2c "bool_pct_overlap" Calculate boolean of overlap % based on given cutoff parameter. More strict than original boolean overlap

# 3 "boolean_distance" ## Calculates distance from edge of assayed region to edge of each chromatin state 
      ###Note, distance measures 3 and 4 are slow, takes 20 minutes for a 200k region file
# 4 "dist_to_feature" , the bp edge to edge data. Requires option 3 calc be ran


#################
### Example usage
# Sample_df <- fread("E2g_file_location.txt")
# source ("projectdir/2025.12_E2G_chromhmm_overlap_function_v0.2.R")
# final.table <- CreateStateIntersectFile(Sample_file= Sample_df,
#                                         ExpName="ATAC_Chmm_Overlap",
#                                         working_directory="/u/projectdir/2025_E2G_ChromHMM_Boolean_Features/",
#                                         DownloadVuModel=F,
#                                         ProvideModel=chromhmm.table, ## or Model_df, Mmust have columns - chr start stop state
#                                         assays_calculating = c("boolean_overlap","bool_pct_overlap","boolean_distance"),
#                                         pct_overlap_cutoff <- 0.5,
#                                         distance_wanted = 10000,
#                                         WriteOutput=T)
#################


CreateStateIntersectFile <- function(Sample_file= Sample_df,
                                     ExpName="ExpName",
                                     working_directory="/projectdir/",
                                     DownloadVuModel=T, ## Download a copr via download.file() to working dir
                                     ProvideModel="none", ## or provide Model_df, Must have columns: chr start end state
                                     assays_calculating = c("boolean_overlap","percent_overlap","majority_state","bool_pct_overlap","boolean_distance","dist_to_feature"),
                                     pct_overlap_cutoff = 0.5, ##Fraction cutoff for "bool_pct_overlap"
                                     distance_wanted = 10000, ##Basepair cutoff for "boolean_distance"
                                     WriteOutput=T){

  print("Starting Calculation")







  
  ##Setup
  ##Folder setup
  outdir.string <- paste0(working_directory,Sys.Date(),"_", ExpName)
  dir.create(outdir.string)
  outdir <- c(paste0(outdir.string,"/"))
  
  

  ### Download are prepare state file
  if(DownloadVuModel==T){
  
    print("Downloading Vu 2022 State Model")
    #### Features to intersect with, prepared with chr start end and some annotation "state"
    ### Vu 100 state chromatin model
    download.file('https://public.hoffman2.idre.ucla.edu/ernst/2K9RS/full_stack/full_stack_annotation_public_release/hg38/hg38_genome_100_segments.bed.gz', destfile = paste0(outdir, "havu_chromhmm_100state.bed.gz"), method = "wget", extra = "-r -p --random-wait")
    chromhmm.table <- fread(paste0(outdir, "havu_chromhmm_100state.bed.gz"), sep="\t")
    colnames(chromhmm.table) <- c("chr", "start", "end", "state")
    chromhmm.table$chromhmm_start <- chromhmm.table$start
    chromhmm.table$chromhmm_end <- chromhmm.table$end
    chromhmm.table$state <- sapply(strsplit(chromhmm.table$state, "_"), "[[",2)
    
    } else {
      chromhmm.table <- ProvideModel
      
      if(ncol(chromhmm.table[,c("chr", "start", "end", "state")])<4){
        print("missing columns (chr/start/end/state) in ProvideModel df, function will fail")
      }
  }
  
  
  
  ######################################################################
  
  
  
  
  ######################################################################
  ###Operation
  
  print("Preparing")
    
  #### Get/reference E2G
  
  e2g.table <- as.data.frame(Sample_file)
  
  ##Prep e2g table for intersection and with an index region name
  tmp.e2g.table <- e2g.table[,c("ElementChr","ElementStart","ElementEnd")]
  tmp.e2g.table$Region <- c(1:nrow(tmp.e2g.table))
  
  tmp.e2g.table$chr <- tmp.e2g.table$ElementChr
  tmp.e2g.table$start <- tmp.e2g.table$ElementStart
  tmp.e2g.table$end <- tmp.e2g.table$ElementEnd
  tmp.e2g.table$peak_start<-tmp.e2g.table$start 
  tmp.e2g.table$peak_end<-tmp.e2g.table$end
  tmp.working.e2g.table <- tmp.e2g.table

  ##Classify states into useful broad categories
  
  
  reference.table <- cbind.data.frame(state=sort(unique(chromhmm.table$state)), Class=NA)
  states.of.interest<-unique(chromhmm.table$state) [grep("EnhA", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Enh_Strong"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("EnhWk", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Enh_Weak"
  
  states.of.interest<- c("EnhA1","EnhA2", "EnhA3", "EnhA4", "EnhA5")
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Enh_Mesenchyme"
  
  states.of.interest<-c("EnhA17","EnhA18", "EnhA19")
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Enh_ESC"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("PromF", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Promoter_Active"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("BivProm", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Promoter_Bivalent"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("Tx", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Transcribed_Region"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("TSS", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "TSS"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("TxEnh", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Enh_Transcribed"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("Acet", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Acetylated_Regions"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("HET", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Heterochromatin_Regions"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("Quies", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "Quiescent_Regions"
  
  states.of.interest<-unique(chromhmm.table$state) [grep("ReprPC", unique(chromhmm.table$state))]
  reference.table[reference.table$state %in% states.of.interest,]$Class <- "RepressivePolycomb_Regions"
  
  
  reference.table$Class[is.na(reference.table$Class)]<- "Other"
  
  
  
  
  
  
  
  
  
  
  ### Build main intersection index
  captured.regions <- genome_intersect(tmp.e2g.table,chromhmm.table, by=c("chr", "start", "end"), mode="both")
  
  
  ## Filter reference to states overlapped
  reference.table <- reference.table[reference.table$state %in% captured.regions$state,]
  
  
  
  
  
  #1 "boolean_overlap" ## Calculates the most simply presence or absence of region to state overlap, given as 1 or 0
  
  ### Initialization statment
  if( length(assays_calculating[assays_calculating %in% c("boolean_overlap","percent_overlap")])>0){
  
    print("Computing boolean_overlap")
    ####Create Boolean Feature columns
  ### Any bp of overlap
  final.column.list <- NULL
  for(Class in unique(reference.table$Class)){
    
    
    states.of.interest <- reference.table[reference.table$Class %in% Class,]$state
    tmp.captured.regions <- captured.regions[captured.regions$state %in% states.of.interest,]
    
    
    tmp.working.e2g.table$Boolean_score <- 0
    tmp.working.e2g.table[tmp.working.e2g.table$Region %in% tmp.captured.regions$Region,]$Boolean_score <- 1
    
    
    colnames(tmp.working.e2g.table)[ncol(tmp.working.e2g.table)] <- paste0("CrhmmBool-", Class )
    
    final.column.list <- c( final.column.list, paste0("CrhmmBool-", Class ))
  }
  
  }
  ## Checks for any overlap
  
  
  #####################################3
  #2 "percent_overlap" ## Calculates the % overlap for each region and state, 0 to 1 (100%), also defines majority state
  
  if( length(assays_calculating[assays_calculating %in% c("percent_overlap","majority_state","bool_pct_overlap")])>0){
    ####Create percent overlap Score for states
    
    final.column.list2<- NULL
    for(Class in unique(reference.table$Class)){
      
      
      states.of.interest <- reference.table[reference.table$Class %in% Class,]$state
      tmp.captured.regions <- captured.regions[captured.regions$state %in% states.of.interest,]
      
      ##Create vectors for comparing overlap
      start.position<- tmp.captured.regions$peak_start
      start.position[(tmp.captured.regions$chromhmm_start - tmp.captured.regions$peak_start)>0] <- tmp.captured.regions$chromhmm_start[(tmp.captured.regions$chromhmm_start - tmp.captured.regions$peak_start)>0]
      
      end.position<- tmp.captured.regions$peak_end
      end.position[(tmp.captured.regions$chromhmm_end - tmp.captured.regions$peak_end)<0] <- tmp.captured.regions$chromhmm_end[(tmp.captured.regions$chromhmm_end - tmp.captured.regions$peak_end)<0]
      
      
      total.peak.size <- tmp.captured.regions$peak_end - tmp.captured.regions$peak_start
      state.overlap.size <- (end.position -start.position)
      
      ##Some states overlap exactly on the start/end bp and give 0% sizes, so adding 1 bp for those
      state.overlap.size[state.overlap.size==0] <-1
      
      
      ##Calculate overlap ratio
      total.state.overlap <- state.overlap.size/total.peak.size
      
      
      ###Need an operator to sum up multi state totals
      ###States are bp specific so can directly sum %
      
      tmp.captured.regions$total_state_overlap <- total.state.overlap
      tmp.captured.regions %>% group_by(Region) %>% summarize(sum=sum(total_state_overlap)) -> overlap.sum.df
      
      tmp.working.e2g.table$Overlap_score <- 0
      tmp.working.e2g.table[tmp.working.e2g.table$Region %in% tmp.captured.regions$Region,]$Overlap_score <- overlap.sum.df$sum
      
      colnames(tmp.working.e2g.table)[ncol(tmp.working.e2g.table)] <- paste0("CrhmmPct-", Class )
      final.column.list2 <- c( final.column.list2, paste0("CrhmmPct-", Class ))
      
    }
    
    ###### 2b Calculate winning max % feature
    ###How many columns in pct.feature?
    total.columns <- length(final.column.list[grep("CrhmmPct-", final.column.list2)] )
    
    tmp.majority.df<-tmp.working.e2g.table[,c( (ncol(tmp.working.e2g.table) - total.columns): ncol(tmp.working.e2g.table))]
    majority.state.vector <- colnames(tmp.majority.df)[apply(tmp.majority.df,1, which.max)]
    majority.state.vector <- gsub("CrhmmPct-","", majority.state.vector)
    
    tmp.working.e2g.table$Majority_State <- majority.state.vector
    
  }
  
  
  ##2c "bool_pct_overlap" Calculate boolean of overlap % based on given cutoff parameter. More strict than original boolean overlap
  
  if( length(assays_calculating[assays_calculating %in% c("bool_pct_overlap")])>0){
    bool_pct_overlap_df<-tmp.working.e2g.table[,final.column.list2]
    bool_pct_overlap_df[bool_pct_overlap_df>pct_overlap_cutoff] = 1
    bool_pct_overlap_df[bool_pct_overlap_df<=pct_overlap_cutoff] = 0
  }
  
  
  
  
  
  ##################################################################
  #3 "boolean_distance" ## Calculates distance from edge of assayed region to edge of each chromatin state 
  
  if( length(assays_calculating[assays_calculating %in% c("boolean_distance","dist_to_feature")])>0){
   
    print("Computing distance matrix")
    
  
  ###############
  ### Module for calculating distance to nearest state
  
  ### Functions for building distance matrix from start and end edges of region to state

  quick_distance_function<-function(x){
    tmp_value <- center.vector[x]
    min(min(abs(tmp_value-tmp_start_vector)) ,min(abs(tmp_value-tmp_end_vector)))
  }
  
  
  ################
  store.df <- NULL
  for(Chromosome in unique(tmp.working.e2g.table$chr)){
    print(Chromosome)
    tmp.chr.states <- chromhmm.table[chromhmm.table$chr==Chromosome,]
    
    tmp.chr.e2g.table <- tmp.working.e2g.table[tmp.working.e2g.table$chr==Chromosome,]
    tmp.chr.e2g.table$center <- (tmp.chr.e2g.table$start+tmp.chr.e2g.table$end)/2
    tmp.chr.e2g.table$size <- abs(tmp.chr.e2g.table$start - tmp.chr.e2g.table$end)
    
    for(Class in unique(reference.table$Class)){
      # print(Class)
      states.of.interest <- reference.table[reference.table$Class %in% Class,]$state
      tmp.chr.states2 <- tmp.chr.states[tmp.chr.states$state %in% states.of.interest,]
      if( nrow(tmp.chr.states2)<1 ){next}
      
      tmp_start_vector<-tmp.chr.states2$start 
      tmp_end_vector<-tmp.chr.states2$end
      center.vector <- tmp.chr.e2g.table$center
      size.vector <- tmp.chr.e2g.table$size/2
      
      tmp.chr.e2g.table$State_Dist<- sapply( c(1:nrow(tmp.chr.e2g.table)),quick_distance_function )-size.vector
      
      tmp.store.df <-tmp.chr.e2g.table[,c("Region","State_Dist")]
      tmp.store.df$Class <- Class
      store.df <- rbind(store.df, tmp.store.df)
      
      
    } 
    
    
  }
    
    
  
  
  
  
  
  
  
  ###Based on this matrix, can make a distance Matrix
  
  tmp.dist.e2g <- dcast(store.df, Region~Class, value.var="State_Dist")
  rownames(tmp.dist.e2g) <- tmp.dist.e2g$Region
  tmp.dist.e2g$Region <- NULL
  
  final.absolute.distance.edge.to.edge <- tmp.dist.e2g
  final.absolute.distance.edge.to.edge[final.absolute.distance.edge.to.edge<0] <- 0
  
  
  
  
  
  
  ###########
  ### Calculate Distance booleans matrix
  
  #### uses dist.wanted parameter
  
  
  bool.dist.region <- -final.absolute.distance.edge.to.edge
  bool.dist.region[bool.dist.region>= (-distance_wanted)] <- 1
  bool.dist.region[bool.dist.region< (-distance_wanted)] <- 0
  
  colnames(bool.dist.region) <- paste0("WithIn",distance_wanted,"_",colnames(bool.dist.region))
  
  
  
  }
  
  
  
  
  
  
  ## If logic to bind based on assays_calculating
  #### rbind just the feature columns, eliminates the pseudo start/end columns I used for calc
  ## 
  
  final.e2g.table <- e2g.table
  
  ##One of first two assays run
  if( length(assays_calculating[assays_calculating %in% c("boolean_overlap","percent_overlap")])>0){
    final.e2g.table <- cbind( final.e2g.table,tmp.working.e2g.table[,final.column.list])
  }
  
  
  ##2b
  if( length(assays_calculating[assays_calculating %in% c("majority_state")])>0){
    final.e2g.table$Majority_State <- tmp.working.e2g.table$Majority_State
  }
  
  ## 2c
  if( length(assays_calculating[assays_calculating %in% c("bool_pct_overlap")])>0){
    final.e2g.table <- cbind(final.e2g.table ,bool_pct_overlap_df)
  }
  
  
  ##3
  if( length(assays_calculating[assays_calculating %in% c("boolean_distance")])>0){
    final.e2g.table<-cbind(final.e2g.table, bool.dist.region)
  }
  
  
  ##4
  if( length(assays_calculating[assays_calculating %in% c("dist_to_feature")])>0){
    tmp.dist.table <- final.absolute.distance.edge.to.edge
    colnames(tmp.dist.table) <- paste0("DistTo_", colnames(tmp.dist.table))
    final.e2g.table<-cbind(final.e2g.table, tmp.dist.table)
  }
  
  
  ###Remove calc columns
  final.e2g.table$peak_start <- NULL
  final.e2g.table$peak_end<- NULL
  final.e2g.table$Region<- NULL
  
    if(WriteOutput==T){
      print("Writing output")
    
    write.table(final.e2g.table, gzfile(paste0(outdir,Exp.Name, "E2G_table.with.ChromHMM.state.overlap.features.txt.gz")),
                sep="\t", quote=F)
    }
  
    final.e2g.table
  }
  




