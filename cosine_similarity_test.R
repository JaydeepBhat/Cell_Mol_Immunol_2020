# Author: Jaydeep Bhat
# Affiliations: 1) Institute of Immunology, Christian-Albrechts-University Kiel & University Hospital Schleswig-Holstein, Campus Kiel, Kiel, Germany; 
# 2) Metabolic Programming, School of Life Sciences Weihenstephan, Technical University Munich (TUM), 85354 Freising, Germany
# Originally developed by Florian Schmidt and Marcel Schulz published in Durek et al (PMID:27851915; see supplemental information). Thanks a lot for sharing the script!
# Acknowledgement: Thanks to Daniela Esser (Institute of Experimental Medicine, Christian-Albrechts-University Kiel & University Hospital Schleswig-Holstein, Campus Kiel, Kiel, Germany) for the help with fixing a minor issue.

# Cosine similarity test for identifying the healthy immune cell signatures in the HSTL patient
HSTL_Healthy <- read.delim("CpGList_HSTL_healthy.txt", header = T, sep = " ")
HSTL_Healthy_Cos0 <- HSTL_Healthy[,c(1:22)] # select sample columns from HSTL patient and healthy T cells
row.names(HSTL_Healthy_Cos0) <- HSTL_Healthy[,1]
numCpGs=dim(HSTL_Healthy_Cos0)[1]
#Initialise lists
values1=c()
values2=c()
values3=c()
values4=c()
values5=c()

#For 100,000 runs
for (i in 1:100000)
{
  
  mod_test <- i%%100
  if(mod_test==0){
    print(paste("done",i,"times",sep=" "))
  }
# sample 266 out of the total numCpGs in the expression data, following the bootstrap, it is with replacement
  All2 < -HSTL_Healthy_Cos0[sample(numCpGs,266,replace=T),]
  All2 <- All2[,2:ncol(All2)]
  All2 <- apply(All2,2,as.numeric)
# compute the mean expression over the replicates
  vis1 <- All2[,1]
  vis2 <- All2[,2]
  vis3 <- All2[,3]
  vis4 <- All2[,4]
  vis5 <- All2[,5]
  Tcells <- apply(All2[,6:21],1,mean)
# compute cosine similiarity score
# single
  Vis1_Tcells <-simil(list(vis1,Tcells),method="cosine")
  Vis2_Tcells <-simil(list(vis2,Tcells),method="cosine")
  Vis3_Tcells <-simil(list(vis3,Tcells),method="cosine")
  Vis4_Tcells <-simil(list(vis4,Tcells),method="cosine")
  Vis5_Tcells <-simil(list(vis5,Tcells),method="cosine")

  values1 <-c(values11,Vis1_Tcells)
  values2 <-c(values12,Vis2_Tcells)
  values3 <-c(values13,Vis3_Tcells)
  values4 <-c(values14,Vis4_Tcells)
  values5 <-c(values15,Vis5_Tcells)
}

# Generate a figure
# All HSTL patient visits and healthy Tcells
boxplotData_v2<-cbind(values1,values2,values3,values4,values5)
save(boxplotData_v2,file="cosine_similarity_total_comparison_3.R")
boxplotData_v2_2<-boxplotData_v2
par(mar=c(12,5,2,2))
pdf("sim_Vis_Tcells_combinations.pdf")
sim_Vis_Tcells_combinations<-boxplot(boxplotData_v2_2,notch=TRUE,las=2,
                                     names=c("Vis1_Tcells","Vis2_Tcells","Vis3_Tcells","Vis4_Tcells","Vis5_Tcells"),
                                     ylim=c(0.91,0.99), ylab="Similiarity score")
boxplotData_order <- c("Vis1_Tcells","Vis2_Tcells","Vis3_Tcells","Vis4_Tcells","Vis5_Tcells")
print(sim_Vis_Tcells_combinations)
dev.off()

#Statistical significance
adjpvalVis12T <- c()
Sim_Vis12_Tcells <- wilcox.test(values1,values2, paired = T)
adjpvalVis12T <- c(adjpvalVis12T,Sim_Vis12_Tcells$p.value)
adjpvalVis12T <- p.adjust(adjpvalVis12T,method = 'fdr')
capture.output(Sim_Vis12_Tcells, file = "Sim_Vis12_Tcells.txt")
Sim_Vis13_Tcells <- wilcox.test(values1,values3, paired = T)
Sim_Vis13_Tcells
capture.output(Sim_Vis13_Tcells, file = "Sim_Vis13_Tcells.txt")
Sim_Vis14_Tcells <- wilcox.test(values1,values4, paired = T)
Sim_Vis14_Tcells
capture.output(Sim_Vis14_Tcells, file = "Sim_Vis14_Tcells.txt")
Sim_Vis15_Tcells <- wilcox.test(values1,values5, paired = T)
Sim_Vis15_Tcells
capture.output(Sim_Vis15_Tcells, file = "Sim_Vis15_Tcells.txt")

sessionInfo()