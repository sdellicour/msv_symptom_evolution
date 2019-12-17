library(ape)
library(ade4)
library(coda)
library(doMC)
library(gplots)
library(raster)
library(vioplot)
library(adephylo)
library(phytools)
library(RColorBrewer)

source("readNexus.r")

logit = function(v)
	{
		for (i in 1:length(v)) v[i] = log10(v[i]/(1-v[i]))
		return(v)
	}
delogit = function(v)
	{
		for (i in 1:length(v)) v[i] = (10^v[i])/(1+(10^v[i]))
		return(v)
	}

wd1 = getwd(); wd2 = paste0(wd1,"/BEAST runs")


# 1. Preparation of the different input files

	# 1.1. Preparation of sequence files

setwd(wd2)
fasta = read.dna("Minus_recombinant_regions_sequences.fasta", format="fasta")
row.names(fasta) = gsub("_","-",row.names(fasta))
row.names(fasta) = gsub("-19","_19",row.names(fasta))
row.names(fasta) = gsub("-20","_20",row.names(fasta))
row.names(fasta) = gsub("-MRCA","_MRCA",row.names(fasta))
for (i in 1:length(row.names(fasta)))
	{
		suffix = unlist(strsplit(row.names(fasta)[i],"_"))[length(unlist(strsplit(row.names(fasta)[i],"_")))]
		if ((suffix == "MRCA")|(suffix == "MRCA-1822-2697"))
			{
				row.names(fasta)[i] = paste("zzzz",row.names(fasta)[i],sep="-")
			}
	}
fasta = fasta[order(row.names(fasta)),]
sequence_names = row.names(fasta)
	
	# 1.2. Estimation and preparation of average values for each individual, trait, host and leaf

		# 1.2.1. Analysis of the minimum and maximal values for each trait
				
setwd(wd1)
tab = read.table("MSV_traits_1.txt", header=T)
for (j in 3:dim(tab)[2])
	{
		t = as.character(tab[1,j])
		t = gsub("c\\(","",t)
		t = gsub("\\)","",t)
		t = unlist(strsplit(t,","))
		t = as.numeric(t)
		minV = min(t); maxV = max(t)
		for (i in 2:dim(tab)[1])
			{
				t = as.character(tab[i,j])
				t = gsub("c\\(","",t)
				t = gsub("\\)","",t)
				t = unlist(strsplit(t,","))
				t = as.numeric(t)
				if (length(t) > 0)
					{
						if (min(t, na.rm=T) < minV)
							{
								minV = min(t, na.rm=T)
							}	
						if (max(t, na.rm=T) > maxV)
							{
								maxV = max(t, na.rm=T)
							}
					}
			}
		cat(paste0("[",round(minV,3),",",round(maxV,3),"]\t\t",colnames(tab)[j]))
		cat("\n")
	}

		# 1.2.2. Removing outliers (Tukey's method) + logit or log_e transformation + standardisation (= always the last step !!)

				# N.B.: boxplot.stats()$out uses the Tukey’s method to identify the outliers ranged above and below the 1.5*IQR
				
logitTransformations = c(F,T,F,T,T,T)
       rescalingStep = c(F,T,F,T,F,F)
minMaxValues = list()
minMaxValues[[2]] = c(0,100)
minMaxValues[[4]] = c(0,50)
traits = list()
averages = list()
logitAverages = list()
CVs = list()
logitCVs = list()
nberOfMeasures = list()
new_ts = tab
new_ts[,3:dim(new_ts)[2]] = -9999
logits = tab
logits[,3:dim(logits)[2]] = -9999
for (i in 3:dim(tab)[2])
	{
		trait = list()
		average = list()
		logitAverage = list()
		CV = list()
		logitCV = list()
		nberOfMeasure = list()
		vS = c()
		for (j in 1:dim(tab)[1])
			{
				t = as.character(tab[j,i])
				t = gsub("c\\(","",t)
				t = gsub("\\)","",t)
				t = unlist(strsplit(t,","))
				t = as.numeric(t)
				vS = c(vS, t)
				outs = boxplot.stats(t)$out
				if (length(outs) > 0)
					{
						t = t[!t[]%in%outs]
					}
				nberOfMeasure[[j]] = length(t)
				trait[[j]] = t
				m = mean(t[!is.na(t[])])
				average[[j]] = m
				CV[[j]] = sd(t[!is.na(t[])])/m
				if (length(t[!is.na(t[])]) == 0)
					{
						text = "c\\(\\)"
					}	else		{
						text = paste("c\\(",t[1],sep="")
						for (k in 2:length(t))
							{
								text = paste(text,t[k],sep=",")
							}
						text = paste(text,"\\)",sep="")
					}
				new_ts[j,i] = text
				if (logitTransformations[i-2] == TRUE)
					{
						if (rescalingStep[i-2] == TRUE)
							{
								minV = minMaxValues[[i-2]][1]
								maxV = minMaxValues[[i-2]][2]
								t = (t-minV)/(maxV-minV) 
							}
						t = logit(t)
					}	else		{
						t = log(t)
					}
				m = mean(t[!is.na(t[])])
				logitAverage[[j]] = m
				logitCV[[j]] = sd(t[!is.na(t[])])/m
				if (length(t[!is.na(t[])]) == 0)
					{
						text = "c\\(\\)"
					}	else		{
						text = paste("c\\(",t[1],sep="")
						for (k in 2:length(t))
							{
								text = paste(text,t[k],sep=",")
							}
						text = paste(text,"\\)",sep="")
					}
				logits[j,i] = text
			}
		traits[[i-2]] = trait
		averages[[i-2]] = average
		logitAverages[[i-2]] = logitAverage
		CVs[[i-2]] = CV
		logitCVs[[i-2]] = logitCV
		nberOfMeasures[[i-2]] = nberOfMeasure
	}
write.table(new_ts, "MSV_traits_2.txt", col.names=T, row.names=F, quote=F, sep="\t")
write.table(logits, "MSV_logits_2.txt", col.names=T, row.names=F, quote=F, sep="\t")

setwd(wd1)
# First: removing all the "\" in a text editor
tab = read.table("MSV_logits_2.txt", header=T)
traitMeans = matrix(nrow=6, ncol=1); traitSDs = matrix(nrow=6, ncol=1)
for (j in 3:dim(tab)[2])
	{
		t = as.character(tab[1,j])
		t = gsub("c\\(","",t)
		t = gsub("\\)","",t)
		t = unlist(strsplit(t,","))
		t = as.numeric(t)
		minV = min(t); maxV = max(t)
		tS = t
		for (i in 2:dim(tab)[1])
			{
				t = as.character(tab[i,j])
				t = gsub("c\\(","",t)
				t = gsub("\\)","",t)
				t = unlist(strsplit(t,","))
				t = as.numeric(t)
				if (length(t) > 0)
					{
						if (min(t, na.rm=T) < minV)
							{
								minV = min(t, na.rm=T)
							}	
						if (max(t, na.rm=T) > maxV)
							{
								maxV = max(t, na.rm=T)
							}
					}
				tS = c(tS, t)
			}
		traitMeans[j-2,1] = mean(tS, na.rm=T); traitSDs[j-2,1] = sd(tS, na.rm=T)
		print(paste(colnames(tab)[j],round(minV,3),round(maxV,3),round(traitMeans[j-2,1],3),round(traitSDs[j-2,1],3),sep=" "))
	}
standAverages = list()
stands = tab
stands[,3:dim(stands)[2]] = -9999
for (i in 3:dim(tab)[2])
	{
		standAverage = list()
		vS = c()
		for (j in 1:dim(tab)[1])
			{
				t = as.character(tab[j,i])
				t = gsub("c\\(","",t)
				t = gsub("\\)","",t)
				t = unlist(strsplit(t,","))
				t = as.numeric(t) # standardisation:
				t = (t-traitMeans[i-2,1])/traitSDs[i-2,1]
				m = mean(t[!is.na(t[])])
				standAverage[[j]] = m
				if (length(t[!is.na(t[])]) == 0)
					{
						text = "c\\(\\)"
					}	else		{
						text = paste("c\\(",t[1],sep="")
						for (k in 2:length(t))
							{
								text = paste(text,t[k],sep=",")
							}
						text = paste(text,"\\)",sep="")
					}
				stands[j,i] = text
			}
		standAverages[[i-2]] = standAverage
	}
write.table(stands, "MSV_stands_2.txt", col.names=T, row.names=F, quote=F, sep="\t")

individuals = unique(tab[,1])
averaged_values_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
averaged_values_2 = matrix(nrow=length(individuals), ncol=(6*3))
logitAveraged_values_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
logitAveraged_values_2 = matrix(nrow=length(individuals), ncol=(6*3))
standAveraged_values_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
standAveraged_values_2 = matrix(nrow=length(individuals), ncol=(6*3))
CV_values_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
CV_values_2 = matrix(nrow=length(individuals), ncol=(6*3))
logitCV_values_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
nberOfMeasures_1 = matrix(nrow=length(individuals), ncol=(6*3*3))
nberOfMeasures_2 = matrix(nrow=length(individuals), ncol=(6*3))
measures = colnames(tab)[3:length(tab)]
host_leafs = as.character(unique(tab[,2]))
hosts = c("SC","GB","Pan77")
colNames_1 = c()
colNames_2 = c()
for (i in 1:length(measures))
	{
		for (j in 1:length(host_leafs))
			{
				colNames_1 = c(colNames_1, paste(measures[i],host_leafs[j],sep="_"))
			}
		for (j in 1:length(hosts))
			{
				colNames_2 = c(colNames_2, paste(measures[i],hosts[j],sep="_"))
			}
	}
colnames(averaged_values_1) = colNames_1
colnames(averaged_values_2) = colNames_2
colnames(logitAveraged_values_1) = colNames_1
colnames(logitAveraged_values_2) = colNames_2
colnames(standAveraged_values_1) = colNames_1
colnames(standAveraged_values_2) = colNames_2
colnames(CV_values_1) = colNames_1
colnames(logitCV_values_1) = colNames_1
colnames(nberOfMeasures_1) = colNames_1
colnames(nberOfMeasures_2) = colNames_2
for (i in 1:length(individuals))
	{
		for (j in 1:length(averages))
			{
				for (k in 1:length(host_leafs))
					{
						index = which((tab[,1]==individuals[i])&(tab[,2]==host_leafs[k]))
						col = ((j-1)*9)+k
						averaged_values_1[i,col] = averages[[j]][[index]]
						logitAveraged_values_1[i,col] = logitAverages[[j]][[index]]
						standAveraged_values_1[i,col] = standAverages[[j]][[index]]
						CV_values_1[i,col] = CVs[[j]][[index]]
						logitCV_values_1[i,col] = logitCVs[[j]][[index]]
						nberOfMeasures_1[i,col] = nberOfMeasures[[j]][[index]]
					}
			}
	}
species_names = c()
for (i in 1:length(tab[,2]))
	{
		species_names = c(species_names, unlist(strsplit(as.character(tab[i,2]),"_"))[1])
	}
for (i in 1:length(individuals))
	{
		for (j in 1:length(traits))
			{
				for (k in 1:length(hosts))
					{
						index = which((tab[,1]==individuals[i])&(species_names==hosts[k]))
						values = c(traits[[j]][[index[1]]], traits[[j]][[index[2]]], traits[[j]][[index[3]]])
						means = c(mean(traits[[j]][[index[1]]]), mean(traits[[j]][[index[2]]]), mean(traits[[j]][[index[3]]]))
						sumOfMeasures = nberOfMeasures[[j]][[index[1]]]+nberOfMeasures[[j]][[index[2]]]+nberOfMeasures[[j]][[index[3]]]
						col = ((j-1)*3)+k
						if (logitTransformations[j] == TRUE)
							{
								if (rescalingStep[j] == TRUE)
									{
										minV = minMaxValues[[j]][1]
										maxV = minMaxValues[[j]][2]
										rescaledValues = (values-minV)/(maxV-minV) 
									}	else	{
										rescaledValues = values
									}
								logitValues = logit(rescaledValues)
							}	else	{
								logitValues = log(values)
							}
						standValues = (logitValues-traitMeans[j,1])/traitSDs[j,1]
						averaged_values_2[i,col] = mean(values[!is.na(values)])
						logitAveraged_values_2[i,col] = mean(logitValues[!is.na(logitValues)])
						standAveraged_values_2[i,col] = mean(standValues[!is.na(standValues)])
						nberOfMeasures_2[i,col] = sumOfMeasures
					}			
			}
	}	
individuals = as.character(individuals)
individuals = gsub("_","-",individuals)
individuals = gsub("-19","_19",individuals)
individuals = gsub("-20","_20",individuals)
individuals = gsub("-MRCA","_MRCA",individuals)
for (i in 1:length(individuals))
	{
		suffix = unlist(strsplit(individuals[i],"_"))[length(unlist(strsplit(individuals[i],"_")))]
		if (suffix == "MRCA")
			{
				individuals[i] = paste("zzzz",individuals[i],sep="-")
			}
	}
row.names(averaged_values_1) = individuals
row.names(averaged_values_2) = individuals
averaged_values_1 = averaged_values_1[order(row.names(averaged_values_1)),]
averaged_values_2 = averaged_values_2[order(row.names(averaged_values_2)),]
row.names(logitAveraged_values_1) = individuals
row.names(logitAveraged_values_2) = individuals
logitAveraged_values_1 = logitAveraged_values_1[order(row.names(logitAveraged_values_1)),]
logitAveraged_values_2 = logitAveraged_values_2[order(row.names(logitAveraged_values_2)),]
row.names(standAveraged_values_1) = individuals
row.names(standAveraged_values_2) = individuals
standAveraged_values_1 = standAveraged_values_1[order(row.names(standAveraged_values_1)),]
standAveraged_values_2 = standAveraged_values_2[order(row.names(standAveraged_values_2)),]
row.names(CV_values_1) = individuals
CV_values_1 = CV_values_1[order(row.names(CV_values_1)),]
row.names(logitCV_values_1) = individuals
logitCV_values_1 = logitCV_values_1[order(row.names(logitCV_values_1)),]
row.names(nberOfMeasures_1) = individuals
row.names(nberOfMeasures_2) = individuals
nberOfMeasures_1 = nberOfMeasures_1[order(row.names(nberOfMeasures_1)),]
nberOfMeasures_2 = nberOfMeasures_2[order(row.names(nberOfMeasures_2)),]

	# 1.3. Saving all the files (mirror versions)

		# 1.3.1. To keep only the first part of each sequence name (removing the sampling year and recombination information)

for (i in 1:length(sequence_names))
	{
		for (j in 1:length(sequence_names[[i]]))
			{
				sequence_names[[i]][j] = unlist(strsplit(sequence_names[[i]][j],"_"))[1]
			}
	}
measure_names = c()
for (i in 1:length(row.names(averaged_values_1)))
	{
		measure_names = c(measure_names, unlist(strsplit(row.names(averaged_values_1)[i],"_"))[1])
	}

		# 1.3.2. To create mirror measure files without individuals with no sequences

trait_indices = c()
nber_of_MRCA = 0
if (length(row.names(fasta)) > 0)
	{
		for (j in 1:length(row.names(fasta)))
			{
				fasta_name = unlist(strsplit(row.names(fasta)[j],"_"))[1]
				trait_indices = c(trait_indices, which(measure_names==fasta_name))
				if (length(grep("zzzz",row.names(fasta)[j])) == 1) nber_of_MRCA = nber_of_MRCA + 1
			}
	}
if (length(names(fasta)) > 0)
	{
		for (j in 1:length(names(fasta)))
			{
				fasta_name = unlist(strsplit(names(fasta)[j],"_"))[1]
				trait_indices = c(trait_indices, which(measure_names==fasta_name))
				if (length(grep("zzzz",names(fasta)[j])) == 1) nber_of_MRCA = nber_of_MRCA + 1
			}
	}
temp = averaged_values_1[trait_indices,]
write.table(temp, "Minus_recombinant_regions_meanValues_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_logitMeanValues_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_standMeanValues_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_CVValues_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_logitCVValues_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_nberOfMeasures_1.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_meanValues_2.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_logitMeanValues_2.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_standMeanValues_2.txt", quote=F)
write.table(temp, "Minus_recombinant_regions_nberOfMeasures_2.txt", quote=F)


# 2. Analysing the correlation among traits

setwd(wd2)
log = read.table("Minus_recombinant_regions_constraints_3.log", header=T)
burnIn = 75000000; log = log[(which(log[,"state"]==burnIn)+1):dim(log)[1],]
hosts = c("SC","GB","Pan77"); host_names = c("M","S","R"); cols = c("#d95f02","#7570b3","#1b9e77")
measures = c("stunting","chlorotic_area","av_streak_width","deviation_from_rectangle","yellowness","whitness")
measure_names = c("leaf stunting","chlorotic area","average streak width","leaf deformation","leaf yellowness","intensity of chlorosis")
indices = which(grepl("correlation",colnames(log))); c = 0
mat1a = matrix(nrow=15, ncol=5); colnames(mat1a) = c(hosts, "trait1", "trait2")
mat1b = matrix(nrow=15, ncol=5); colnames(mat1b) = c(hosts, "trait1", "trait2")
for (i in 1:length(hosts))
	{
		l = 0
		for (j in 1:(length(measures)-1))
			{
				for (k in (j+1):length(measures))
					{
						l = l+1; c = c+1; mat1a[l,i] = mean(log[,indices[c]])
						quantiles = quantile(log[,indices[c]], c(0.025,0.975))
						significant = TRUE
						if ((quantiles[1]<0)&(0<quantiles[2])) significant = FALSE
						if (significant == FALSE) mat1b[l,i] = 0
						if (significant == TRUE) mat1b[l,i] = 1
						if (i == 1)
							{
								mat1a[l,"trait1"] = measures[j]; mat1b[l,"trait1"] = measures[j]
								mat1a[l,"trait2"] = measures[k]; mat1b[l,"trait2"] = measures[k]
							}
					}
			}
	}
hosts = c("Pan77","SC","GB"); host_names = c("R","M","S"); cols = c("#1b9e77","#d95f02","#7570b3")
measures = c("stunting","chlorotic_area","deviation_from_rectangle","whitness")
measure_names = c("leaf stunting","chlorotic area","leaf deformation","intensity of chlorosis")
mat2a = matrix(nrow=6, ncol=3); colnames(mat2a) = hosts
mat2b = matrix(nrow=6, ncol=3); colnames(mat2b) = hosts
rowNames = rep(NA, dim(mat2a)[1])
for (i in 1:length(hosts))
	{
		l = 0
		for (j in 1:(length(measures)-1))
			{
				for (k in (j+1):length(measures))
					{
						index1 = which((mat1a[,"trait1"]==measures[j])&(mat1a[,"trait2"]==measures[k]))
						index2 = which(colnames(mat1a)==hosts[i]); l = l+1
						mat2a[l,i] = as.numeric(mat1a[index1,index2]); mat2b[l,i] = as.numeric(mat1b[index1,index2])
						if (i == 1) rowNames[l] = paste0(measure_names[j]," - ",measure_names[k])
					}
			}
	}
row.names(mat2a) = rowNames; row.names(mat2b) = rowNames; mat2a[1:3,] = -mat2a[1:3,]
dev.new(width=5, height=3.5)
heatmap.2(mat2a, main="", density.info="none", dendrogram="none", srtCol=45, key=T, cellnote=mat2c, notecex=0.6, notecol="black",
labRow=row.names(mat2a), labCol=host_names, cexRow=0.8, cexCol=0.9, cex=1, sepcolor="black", colsep=c(0:10), rowsep=c(0:15),
sepwidth=c(0.001, 0.001), margin=c(8,15), trace="none", col=colorRampPalette(brewer.pal(11,"RdYlGn"))(101), Rowv=NULL, Colv="NA")
print(mat2a); print(mat2b); # dev.off()


# 3. Analysing and comparing phylogenetic signals

setwd(wd2)
log = read.table("Minus_recombinant_regions_constraints_2.log", header=T)
burnIn = 101; log = log[(burnIn+1):dim(log)[1],]
measures = c("av_streak_width","chlorotic_area","deviation_from_rectangle","stunting","yellowness","whitness")
hosts = c("Pan77","SC","GB")
meanLambdas = matrix(nrow=6, ncol=3); colnames(meanLambdas) = hosts; rownames(meanLambdas) = measures
hpdLambdas = matrix(nrow=6, ncol=3); colnames(hpdLambdas) = hosts; rownames(hpdLambdas) = measures
for (i in 1:length(measures))
	{
		for (j in 1:length(hosts))
				{
					col = which(grepl("lambda", colnames(log))==TRUE)
					
					lambda = log[,paste(measures[i],"_",hosts[j],".lambda",sep="")]
					meanLambda = mean(lambda)
					hpdLambda = HPDinterval(as.mcmc(lambda), 0.95)
					hpdLambda = paste(round(hpdLambda[1,1],3),"-",round(hpdLambda[1,2],3),sep="")
					meanLambdas[i,j] = round(meanLambda,3)
					hpdLambdas[i,j] = hpdLambda
			}
	}
write.table(meanLambdas, "Minus_recombinant_regions_mean_lambdas.txt", quote=F, sep="\t")
write.table(hpdLambdas, "Minus_recombinant_regions_HPD_lambdas.txt", quote=F, sep="\t")


# 4. Continuous character mapping analyses (i.e. plotting coloured MCC trees)

	# 4.1. Reading and transforming back all the inferred values

setwd(wd2)
meanLambdas = read.table("Minus_recombinant_regions_mean_lambdas.txt", header=T)
hpdLambdas = read.table("Minus_recombinant_regions_HPD_lambdas.txt", header=T)
tree1 = readNexus("Minus_recombinant_regions_constraints_MCC1.tree"); tree2 = tree1
nodesToRotate = c(119,97,102,106,110,103,105,111,117,107,64,77,78,82,84,86,92,95,88,89,90,85,75,67,69)
for (i in nodesToRotate) tree2 = ape::rotate(tree2, node=i) # CAUTION: the annotations are not rotated!
annotations = list()
indicesAfterRotations = matrix(nrow=dim(tree2$edge)[1], ncol=1)
for (i in 1:dim(indicesAfterRotations)[1])
	{
		indicesAfterRotations[i,1] = which((tree1$edge[,1]==tree2$edge[i,1])&(tree1$edge[,2]==tree2$edge[i,2]))
	}
tree = tree2
for (i in 1:length(indicesAfterRotations))
	{
		annotations[[i]] = tree1$annotations[[indicesAfterRotations[i]]]
	}
originalTipLabels = tree$tip.label; nberOfTipNodes = length(originalTipLabels)
indexOfTipNodes = matrix(nrow=nberOfTipNodes, ncol=1)
for (i in 1:nberOfTipNodes)
	{
		indexOfTipNodes[i,1] = which(tree$edge[,2]==i)
	}
indexOfInternalNodes = matrix(nrow=tree$Nnode, ncol=1)
for (i in 1:tree$Nnode)
	{
		if (i > 1) indexOfInternalNodes[i,1] = which(tree$edge[,2]==(nberOfTipNodes+i))
	}
measures = c("av_streak_width","deviation_from_rectangle","chlorotic_area","stunting","whitness","yellowness")
logitTransformations = c(F,T,T,F,T,T)
       rescalingStep = c(F,T,T,F,F,F)      
	   minMaxValues = list()
	   minMaxValues[[3]] = c(0,100)
	   minMaxValues[[2]] = c(0,50)
hosts = c("GB","SC","Pan77"); colNames = c()
for (i in 1:length(measures))
	{
		for (j in 1:length(hosts)) colNames = c(colNames, paste(measures[i],hosts[j],sep="_"))
	}
means = matrix(nrow=length(annotations), ncol=3*6); colnames(means) = colNames
hpdlow = matrix(nrow=length(annotations), ncol=3*6); colnames(hpdlow) = colNames
hpdhigh = matrix(nrow=length(annotations), ncol=3*6); colnames(hpdhigh) = colNames
rates = matrix(nrow=length(annotations), ncol=3*6); colnames(rates) = colNames 
for (i in 1:dim(means)[1])
	{
		for (j in 1:length(measures))
			{
				for (k in 1:length(hosts))
					{
						trait = paste(measures[j],hosts[k],sep="_")
						index = which(names(annotations[[i]])==trait)
						value1 = as.numeric(annotations[[i]][index])
						if (logitTransformations[j] == TRUE)
							{
								value2 = delogit(value1)
								if (rescalingStep[j] == TRUE)
									{
										minV = minMaxValues[[j]][1]
										maxV = minMaxValues[[j]][2]
										value2 = (value2*(maxV-minV))+minV
									}
							}
						if (logitTransformations[j] == FALSE)
							{
								value2 = exp(value1)
							}
						means[i,((j-1)*3)+k] = value2
						trait = paste(measures[j],hosts[k],"95%_HPD",sep="_")
						if (trait%in%names(annotations[[i]]))
							{
								index = which(names(annotations[[i]])==trait)
								value1_low = as.numeric(annotations[[i]][index][[1]])[1]
								value1_high = as.numeric(annotations[[i]][index][[1]])[2]
								if (logitTransformations[j] == TRUE)
									{
										value2_low = delogit(value1_low)
										value2_high = delogit(value1_high)
										if (rescalingStep[j] == TRUE)
											{
												minV = minMaxValues[[j]][1]
												maxV = minMaxValues[[j]][2]
												value2_low = (value2_low*(maxV-minV))+minV
												value2_high = (value2_high*(maxV-minV))+minV
											}
									}
								if (logitTransformations[j] == FALSE)
									{
										value2_low = exp(value1_low)
										value2_high = exp(value1_high)
									}
								hpdlow[i,((j-1)*3)+k] = value2_low
								hpdhigh[i,((j-1)*3)+k] = value2_high
							}
						rate = paste(measures[j],"_",hosts[k],".rate",sep="")
						if (rate%in%names(annotations[[i]]))
							{
								index = which(names(annotations[[i]])==rate)
								value = as.numeric(annotations[[i]][index])
								rates[i,((j-1)*3)+k] = value
							}
					}
			}
	}
tipNodesMeasures = means[indexOfTipNodes,]; row.names(tipNodesMeasures) = originalTipLabels
internalNodesMeasures = means[indexOfInternalNodes,]

	# 4.2. Comparing measured and inferred trait measures for ancestors

setwd(wd2)
tab1 = read.table("Minus_recombinant_regions_meanValues_2.txt", header=T)
tab1 = tab1[which(grepl("zzzz",rownames(tab1))),]
rownames(tab1) = gsub("zzzz_","",gsub("_MRCA","",gsub("-","_",rownames(tab1))))
tab1 = tab1[which(!grepl("o65",rownames(tab1))),]
setwd(wd2)
tree$tip.label = originalTipLabels
tips_A3 = c("KE-O57-Kimathi-est_1990","KE-o54b-Lerosho-kenya_1990","MSV-A-KE-MtKA_1997")
A3 = findMRCA(tree, tips_A3, type="node")
tips_A4 = c("MSV-A-ZA-Cat2-D3_2006","MSV-A-ZA-Kom_1989","MSV-A-ZA-Pot1-Riz48_2007","MSV-A-ZA-Pot4-O28_1979",
		 "MSV-A-ZA-Pot8-O33_1979","MSV-A-ZW-Hel2-Bet36_2006","MSV-A-ZW-Mas1-Bet43_2006","MSV-A4-Za-Nat-g195_2006",
		 "MSV-A4-Za-RosE-g131_2006","MSV-ZA-Pot3-O27K_1979","New-ZA-O59-Greytown_1989","New-ZA-O60-Bethlehem_1989","ZA-O82-letsele_1987")
A4 = findMRCA(tree, tips_A4, type="node")
tips_Letsele_makd = c("MSV-A-MZ-Bil6-Bet25_2007","MSV-A-MZ-Chi1-chimoz_2007","MSV-A-MZ-Inh2-Moz1_2007","MSV-A-MZ-Map8-Moz3_2007",
		 "MSV-A-MZ-Map9-Moz4_2007","MSV-A-MZ-Xai1-xaimoz_2007","MSV-A-UG-Kib188_2005","MSV-A-UG-Luw110_2005","MSV-A-ZA-Mak2-M49_1988",
		 "MSV-A-ZA-MakD_1998","MSV-A1-Za-ThoE-g132_2006","New-ZA-O87RC-letsele_1987")
Letsele_makd = findMRCA(tree, tips_Letsele_makd, type="node")
tips_Sag_mic22 = c("MSV-A-KE-Ama_1998","MSV-A-MZ-Pem2-Moz37_2007","MSV-A-MZ-Pem5-Moz41_2007","MSV-A-UG-Bush53_2005","MSV-A-UG-Hoi154_2005",
		 "MSV-A-UG-Kab82_2005","MSV-A-UG-Kap292_2005","MSV-A-UG-Mba41_2005","MSV-A-UG-Mpi-11_2005","MSV-A-UG-Mub94_2005","MSV-A-ZW-Chi-Zim3_2006",
		 "MSV-A-ZW-Har2-Mic22_1987","MSV-A-ZW-Mas2-Mic4_1993","MSV-A-ZW-MatA_1994","MSV-ULuw-107_2005","MSV-UWak-4_2005")
Sag_mic22 = findMRCA(tree, tips_Sag_mic22, type="node")
tips_A1 = c(tips_Sag_mic22,tips_Letsele_makd,"MSV-A-ZW-MatC_1998","New-ZA-O45RC-Ilanga_1988",
		 "MSV-A-UG-Kib179_2005","MSV-UKas-75_2005","MSV-A-UG-Bug245_2005","New-ZA-O83-letsele_1987")
A1 = findMRCA(tree, tips_A1, type="node")
tips_A1234 = c(tips_A1,tips_A3,tips_A4,"MSV-A-NG-Ns_1980","MSV-A4-Za-Omr-g221_2007","MSV-A4-Zw-Nmg-g168_2006",
		 "MSV-A-ZA-Pot6-O29_1979","New-ZA-O65-Bethlehem_1989","MSV-A-ZA-Hei-O9_1979")
A1234 = findMRCA(tree, tips_A1234, type="node")
tips_A12345 = c(tips_A1234,"MSV-A-RE-Pie4-Mic13_1994","MSV-A-RE-Reu2_1997","MSV-A-RE-Jos1-Mic18_1995")
A12345 = findMRCA(tree, tips_A12345, type="node")
ancestors = c(A1234, A1, A3, A4, Letsele_makd, Sag_mic22, A12345)
internalNodes = seq(61,119,1)
internalNodes[which(!seq(61,119,1)%in%ancestors)] = NA
ancestor_names = c("Anc1","Anc2","Anc3","Anc4","Anc5","Anc6","Anc0")
for (i in 1:length(ancestors))
	{
		internalNodes[which(internalNodes==ancestors[i])] = ancestor_names[i]
	}
tree$node.label = internalNodes
indices = ancestors; indices[] = NA
for (i in 1:length(ancestors))
	{
		indices[i] = which(tree$edge[,2]==ancestors[i])
	}
tab2 = tab1; tab2[,] = NA
tab2_low = tab1; tab2_low[,] = NA
tab2_high = tab1; tab2_high[,] = NA
for (i in 1:length(ancestors))
	{
		for (j in 1:length(colnames(tab2)))
			{
				tab2[i,j] = means[indices[i],colnames(tab2)[j]]
				tab2_low[i,j] = hpdlow[indices[i],colnames(tab2_low)[j]]
				tab2_high[i,j] = hpdhigh[indices[i],colnames(tab2_high)[j]]
			}
	}
tab3 = tab1; tab3[] = NA 
tab4 = tab1; tab4[] = NA
significant_digits = c(3,3,3,1,1,1,4,4,4,3,3,3,3,3,3,3,3,3)
for (i in 1:dim(tab3)[2])
	{
		for (j in 1:dim(tab3)[1])
			{
				meanValue = round(tab2[j,i],significant_digits[i])
				lowerValue = round(tab2_low[j,i],significant_digits[i])
				higherValue = round(tab2_high[j,i],significant_digits[i])
				tab4[j,i] = paste0(meanValue," [",lowerValue,"-",higherValue,"]")
			}
	}
write.csv(tab4, "Mean_95%CI_traits_inferred_for_ancestors.csv", quote=F)

	# 4.3. 95% CI intervals related to ancestor measures (CI on distributions)

setwd(wd1) # First: removing all the "\" in TextWrangler !!!!
tab = read.table("MSV_traits_2.txt", header=T)
indices = which(grepl("-MRCA",tab[,1]))
tab1_low = tab1; tab1_low[,] = NA
tab1_high = tab1; tab1_high[,] = NA
for (i in 1:dim(tab1)[1])
	{
		for (j in 1:dim(tab1)[2])
			{
				tS = c()
				host = unlist(strsplit(colnames(tab1)[j],"_"))[length(unlist(strsplit(colnames(tab1)[j],"_")))]
				lines = which(grepl("-MRCA",tab[,1]) & grepl(paste0(row.names(tab1)[i],"_MRCA"),gsub("-","_",tab[,1])) & grepl(host,tab[,2]))
				column = which(colnames(tab) == gsub(paste("_",host,sep=""),"",colnames(tab1)[j]))
				for (k in 1:length(lines))
					{
						t = as.character(tab[lines[k],column])
						t = gsub("c\\(","",t)
						t = gsub("\\)","",t)
						t = unlist(strsplit(t,","))
						t = as.numeric(t)
						tS = c(tS, t)
					}
				tS = tS[!is.na(tS)]
				quantiles = quantile(tS[!is.na(tS[])],probs=c(0.025,0.975))
				tab1_low[i,j] = quantiles[1]
				tab1_high[i,j] = quantiles[2]
			}
	} 
tab5 = tab1; tab5[] = NA
significant_digits = c(3,3,3,1,1,1,4,4,4,3,3,3,3,3,3,3,3,3)
for (i in 1:dim(tab5)[2])
	{
		for (j in 1:dim(tab5)[1])
			{
				meanValue = round(tab1[j,i],significant_digits[i])
				lowerValue = round(tab1_low[j,i],significant_digits[i])
				higherValue = round(tab1_high[j,i],significant_digits[i])
				tab5[j,i] = paste0(meanValue," [",lowerValue,"-",higherValue,"]")
			}
	}
setwd(wd2)
write.csv(tab5, "Mean_95%CI_traits_observed_for_ancestors1.csv", quote=F) # CI on distributions

	# 4.4. 95% CI intervals related to ancestor measures (CI on means)

setwd(wd1) # First: removing all the "\" in TextWrangler !!!!
tab = read.table("MSV_traits_2.txt", header=T)
indices = which(grepl("-MRCA",tab[,1]))
tab1_mean = tab1; tab1_low[,] = NA
tab1_low = tab1; tab1_low[,] = NA
tab1_high = tab1; tab1_high[,] = NA
for (i in 1:dim(tab1)[1])
	{
		for (j in 1:dim(tab1)[2])
			{
				tS = c()
				host = unlist(strsplit(colnames(tab1)[j],"_"))[length(unlist(strsplit(colnames(tab1)[j],"_")))]
				lines = which(grepl("-MRCA",tab[,1]) & grepl(paste0(row.names(tab1)[i],"_MRCA"),gsub("-","_",tab[,1])) & grepl(host,tab[,2]))
				column = which(colnames(tab) == gsub(paste("_",host,sep=""),"",colnames(tab1)[j]))
				for (k in 1:length(lines))
					{
						t = as.character(tab[lines[k],column])
						t = gsub("c\\(","",t)
						t = gsub("\\)","",t)
						t = unlist(strsplit(t,","))
						t = as.numeric(t)
						tS = c(tS, t)
					}
				tS = tS[!is.na(tS)]
				n = length(tS); m = mean(tS)
				sd = sd(tS) # sample standard deviation
				se = sd(tS)/sqrt(n) # standard error estimate
				error = qt(0.975, df=n-1)*se
				tab1_mean[i,j] = m
				tab1_low[i,j] = m-error
				tab1_high[i,j] = m+error
			}
	} 
tab6 = tab1; tab6[] = NA
significant_digits = c(3,3,3,1,1,1,4,4,4,3,3,3,3,3,3,3,3,3)
for (i in 1:dim(tab5)[2])
	{
		for (j in 1:dim(tab5)[1])
			{
				meanValue = round(tab1_mean[j,i],significant_digits[i])
				lowerValue = round(tab1_low[j,i],significant_digits[i])
				higherValue = round(tab1_high[j,i],significant_digits[i])
				tab6[j,i] = paste0(meanValue," [",lowerValue,"-",higherValue,"]")
			}
	}
setwd(wd2)
write.csv(tab6, "Mean_95%CI_traits_observed_for_ancestors2.csv", quote=F) # CI on means

	# 4.5. Generating and saving all the tree plots with mean values

host_names = c("S","M","R")
measure_names = c("average streak width","leaf deformation","chlorotic area","leaf stunting","intensity of chlorosis","leaf yellowness")
for (i in 1:dim(means)[2])
	{
		if (grepl("stunting",colnames(means)[i])) means[,i] = 1-means[,i]
	}
minMaxValues = matrix(nrow=2, ncol=6); colnames(minMaxValues) = measures; row.names(minMaxValues) = c("min","max")
for (i in 1:dim(minMaxValues)[2])
	{
		columns = which(grepl(colnames(minMaxValues)[i], colnames(means)))
		minMaxValues[1,i] = min(means[,columns]); minMaxValues[2,i] = max(means[,columns])
	}
labelsValues = matrix(nrow=5, ncol=6); colnames(minMaxValues) = measures; nberOfDecimals = c(5,2,2,2,2,2)
for (i in 1:dim(minMaxValues)[2])
	{
		columns = which(grepl(colnames(minMaxValues)[i], colnames(means)))
		labelsValues[1,i] = round(min(means[,columns]), nberOfDecimals[i])
		labelsValues[2,i] = round(min(means[,columns])+ ((max(means[,columns])-min(means[,columns]))*0.25), nberOfDecimals[i])
		labelsValues[3,i] = round(min(means[,columns])+ ((max(means[,columns])-min(means[,columns]))*0.50), nberOfDecimals[i])
		labelsValues[4,i] = round(min(means[,columns])+ ((max(means[,columns])-min(means[,columns]))*0.75), nberOfDecimals[i])
		labelsValues[5,i] = round(max(means[,columns]), nberOfDecimals[i])
	}
colorBrewers = list()
colorBrewers[[1]] = rev(colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121])
colorBrewers[[2]] = colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121]
colorBrewers[[3]] = colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121]
colorBrewers[[4]] = colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121]
colorBrewers[[5]] = colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121]
colorBrewers[[6]] = colorRampPalette(brewer.pal(9,"Greys"))(141)[21:121]
indices = c(2,4,3,5)
for (i in 1:2)
	{
		pdf(paste("Minus_recombinant_regions_MCCtrees_",i,"_means.pdf",sep=""), width=8.00, height=9.00)
		par(mfrow=c(3,2), mar=c(1,1,1,1), oma=c(0,0,0,0), lwd=0.2)
		for (j in 1:3)
			{
				host = hosts[j]
				for (k in 1:2)
					{
						if (k == 1) treeDirection = "rightwards"
						if (k == 2) treeDirection = "leftwards"
						M = ((i-1)*2)+k; M = indices[M]; measure = measures[M]
						index = which(grepl(host,colnames(means))&grepl(measure,colnames(means)))
						cols = colorBrewers[[M]][(((means[,index]-minMaxValues[1,M])/(minMaxValues[2,M]-minMaxValues[1,M]))*100)+1]
						plot(tree, edge.color=cols, edge.width=2.2, direction = treeDirection, show.tip.label=F, show.node.label=T)
						title = paste(gsub("_"," ",measure_names[M]),host_names[j],sep=" - ")
						lambda1 = paste("lambda = ",meanLambdas[measure,host],sep="")
						lambda2 = paste("95% HPD = [",as.character(hpdLambdas[measure,host]),"]",sep="")
						if (j%in%c(1,2,3))
							{
								mat = minMaxValues[,M]; rast = raster(matrix(mat))
								if (k == 1)
									{
										text(x=40, y=35, title, cex=0.9)
										text(x=40, y=32, lambda1, cex=0.8)
										text(x=40, y=29, lambda2, cex=0.8)
										plot(rast, legend.only=T, add=T, col=colorBrewers[[M]], legend.width=0.5, legend.shrink=0.3,
										smallplot=c(0.21,0.23,0.1,0.42), legend.args=list(text="", cex=0.8, line=0.5, col="gray30"),
										axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-0.5, line=0, mgp=c(0,0.5,0)), alpha=1)
									}	else		{
										text(x=205-35, y=35, title, cex=0.9)
										text(x=205-35, y=32, lambda1, cex=0.8)
										text(x=205-35, y=29, lambda2, cex=0.8)
										plot(rast, legend.only=T, add=T, col=colorBrewers[[M]], legend.width=0.5, legend.shrink=0.3,
										smallplot=c(0.73,0.75,0.1,0.42), legend.args=list(text="", cex=0.8, line=0.5, col="gray30"),
										axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-0.5, line=0, mgp=c(0,0.5,0)), alpha=1)
									}
							}
						if (j == 3)
							{
								rootTime = 2017-max(node.depth.edgelength(tree))
								axisPhylo(lwd=0, cex.axis=0.6, mgp=c(0,0.2,-0.2), axis.args=list(at=seq(1800,2000,50)), pos=-2.0,
										  backward=F, root.time=rootTime, lwd.tick=0.2, tck=-0.012, side=1, col="gray30", col.tck="gray30")
								axis(lwd=0.2, at=c(rootTime,0), labels=c("",""), tck=0, pos=-2.0, side=1, col="gray30", col.tck="gray30")
							}
					}
			}
		dev.off()
	}

	# 4.6. Generating graph of the average trait evolution (skyline-like with HPD intervals)

setwd(wd2)
trees = scan("Minus_recombinant_regions_constraints_2.trees", what="", sep="\n", quiet=TRUE)
index0 = which(grepl("tree STATE_0",trees))-1
index1 = which(grepl("tree STATE_100400000",trees))
index2 = which(grepl("tree STATE_500000000",trees))
indices = seq(index1, index2, 4)
trees = c(trees[c(1:index0,indices)],"End;")
write(trees, "Minus_recombinant_regions_constraints_2_1000.trees")

setwd(wd2)
trees = readNexus("Minus_recombinant_regions_constraints_2_1000.trees")

mostRecentSamplingYear = 2007
colNames = c("node1","node2","startYear","endYear"); hosts = c("Pan77","SC","GB")
measures = c("av_streak_width","chlorotic_area","deviation_from_rectangle","stunting","yellowness","whitness")
logitTransformations = c(F,T,T,F,T,T)
       rescalingStep = c(F,T,T,F,F,F)      
	   minMaxValues = list()
	   minMaxValues[[2]] = c(0,100)
	   minMaxValues[[3]] = c(0,50)
for (i in 1:length(hosts))
	{
		for (j in 1:length(measures))
			{
				colNames = c(colNames, paste("start",measures[j],hosts[i],sep="_"), paste("end",measures[j],hosts[i],sep="_"))
			}
	}
for (i in 1:length(trees))
	{
		tree = trees[[i]]
		tab1 = matrix(nrow=dim(tree$edge)[1], ncol=length(colNames))
		tab2 = matrix(nrow=dim(tree$edge)[1], ncol=length(colNames))
		colnames(tab1) = colNames; colnames(tab2) = colNames
		tab1[,c("node1","node2")] = tree$edge; tab2[,c("node1","node2")] = tree$edge
		tab1[,c("startYear","endYear")] = 2007-round((max(nodeHeights(tree))-nodeHeights(tree)),5)
		tab2[,c("startYear","endYear")] = 2007-round((max(nodeHeights(tree))-nodeHeights(tree)),5)
		for (j in 1:dim(tab1)[1])
			{
				annotations = tree$annotations[[j]]
				for (k in 5:dim(tab1)[2])
					{
						if (grepl("end_",colnames(tab1)[k])) tab1[j,k] = as.numeric(annotations[gsub("end_","",colnames(tab1)[k])])
					}
				index = which(tab1[,"node2"]==tab1[j,"node1"])
				if (length(index) > 0)
					{
						annotations = tree$annotations[[index]]
						for (k in 5:dim(tab1)[2])
							{
								if (grepl("start_",colnames(tab1)[k])) tab1[j,k] = as.numeric(annotations[gsub("start_","",colnames(tab1)[k])])
							}
					}
				if (!tab1[j,"node1"]%in%tab1[,"node2"])
					{
						annotations = tree$root.annotation
						for (k in 5:dim(tab1)[2])
							{
								if (grepl("start_",colnames(tab1)[k])) tab1[j,k] = as.numeric(annotations[gsub("start_","",colnames(tab1)[k])])
							}
					}
			}
		ancestralBranches = which(!tab1[,"node1"]%in%tab1[,"node2"])
		tab1 = tab1[-ancestralBranches,]; tab2 = tab2[-ancestralBranches,] # to remove the outgroup ("VW")
		ancestralBranches = which(!tab1[,"node1"]%in%tab1[,"node2"])
		for (j in 5:dim(tab1)[2])
			{
				measure = gsub("start_","",gsub("end_","",colnames(tab1)[j]))
				measure = gsub("_Pan77","",gsub("_SC","",gsub("_GB","",measure)))
				index = which(measures==measure)
				values1 = tab1[,j]
				if (logitTransformations[index] == TRUE)
					{
						values2 = delogit(values1)
						if (rescalingStep[index] == TRUE)
							{
								minV = minMaxValues[[index]][1]
								maxV = minMaxValues[[index]][2]
								values2 = (values2*(maxV-minV))+minV
							}
					}
				if (logitTransformations[index] == FALSE)
					{
						values2 = exp(values1)
					}
				tab2[,j] = values2 
			}
		write.csv(tab2, paste0("Minus_recombinant_regions_constraints_2_1000/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}

traitEvolution = function(localTreesDirectory="", nberOfExtractionFiles=1, timeSlices=100, slidingWindow=1, trait="")
	{
		for (t in 1:nberOfExtractionFiles)
			{
				fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
				data = read.csv(fileName, h=T)
				data = data[with(data, order(endYear,startYear)),]
				nberOfConnections = dim(data)[1]
				if (t == 1)
					{
						minStartYear = min(data[,"startYear"])
			    		minEndYear = min(data[,"startYear"])
			    		maxEndYear = max(data[,"endYear"])
					}	else	{
						if (minStartYear > min(data[,"startYear"])) minStartYear = min(data[,"startYear"])
						if (minEndYear > min(data[,"endYear"])) minEndYear = min(data[,"endYear"])
						if (maxEndYear < max(data[,"endYear"])) maxEndYear = max(data[,"endYear"])
					}
			}
		traitValuesTimeInterval = (maxEndYear-minStartYear)/timeSlices
		xLim = c(minStartYear, maxEndYear)
		traitValuesList = list()
		startEndTimes = matrix(nrow=timeSlices, ncol=3)
		for (i in 1:timeSlices)
			{
				time = minStartYear+((i-1)*traitValuesTimeInterval)+(traitValuesTimeInterval/2)
				startTime = time - (slidingWindow/2)
				endTime = time + (slidingWindow/2)
				startEndTimes[i,1:3] = cbind(time, startTime, endTime)
			}
		buffer = list()
		for (t in 1:nberOfExtractionFiles)
			{
				fileName = paste(localTreesDirectory,"/TreeExtractions_",t,".csv",sep="")
				data = read.csv(fileName, h=T)
				data = data[with(data, order(endYear, startYear)),]
				nberOfConnections = dim(data)[1]
				values = matrix(nrow=timeSlices, ncol=2)
				data = data[order(data[,"endYear"]),]
				for (i in 1:timeSlices)
					{
						vS = c()
						time = startEndTimes[i,1]
						startTime = startEndTimes[i,2]
						endTime = startEndTimes[i,3]
						for (j in 1:nberOfConnections)
							{
								branchInTimeInterval = FALSE
								if ((data[j,"startYear"]<startTime)&(data[j,"endYear"]>endTime)) branchInTimeInterval = TRUE
								if ((data[j,"startYear"]>startTime)&(data[j,"startYear"]<endTime)) branchInTimeInterval = TRUE
								if ((data[j,"endYear"]>startTime)&(data[j,"endYear"]<endTime)) branchInTimeInterval = TRUE
								if (branchInTimeInterval == TRUE)
									{
										if (data[j,"startYear"] > startTime) time1 = data[j,"startYear"]
										if (data[j,"startYear"] <= startTime) time1 = startTime
										if (data[j,"endYear"] < endTime) time2 = data[j,"endYear"]
										if (data[j,"endYear"] >= endTime) time2 = endTime
										timeProportion = (time-data[j,"startYear"])/(data[j,"endYear"]-data[j,"startYear"])
										v = data[j,paste0("start_",trait)]+((data[j,paste0("end_",trait)]-data[j,paste0("start_",trait)])*timeProportion)
			    						vS = c(vS, v)
									}
							}
						if (time > 0)
							{
								values[i,1] = time; values[i,2] = mean(vS)
							}
					}
				colnames(values) = c("year", "traitValues")
				buffer[[t]] = values
			}
		for (t in 1:length(buffer))
			{
				traitValuesList[[t]] = buffer[[t]]
			}
		slicedTimes = matrix(nrow=1, ncol=timeSlices)
		lower_l_95 = matrix(nrow=1, ncol=timeSlices)
		upper_l_95 = matrix(nrow=1, ncol=timeSlices)
		lower_l_80 = matrix(nrow=1, ncol=timeSlices)
		upper_l_80 = matrix(nrow=1, ncol=timeSlices)
		traitValues = matrix(nrow=nberOfExtractionFiles, ncol=timeSlices)
		traitMeanValue = matrix(nrow=1, ncol=timeSlices)
		traitMedianValue = matrix(nrow=1, ncol=timeSlices)
		traitValues = matrix(nrow=nberOfExtractionFiles, ncol=timeSlices)
		traitMeanValue = matrix(nrow=1, ncol=timeSlices)
		traitMedianValue = matrix(nrow=1 ,ncol=timeSlices)
		for (i in 1:timeSlices)
			{
				slicedTimes[1,i] = traitValuesList[[1]][i,1]
				for (t in 1:nberOfExtractionFiles)
					{
						traitValues[t,i] = traitValuesList[[t]][i,2]
					}
				quantiles = quantile(traitValues[,i], probs=c(0.025,0.975), na.rm=T)
				lower_l_95[1,i] = as.numeric(quantiles[1])
				upper_l_95[1,i] = as.numeric(quantiles[2])
				quantiles = quantile(traitValues[,i], probs=c(0.1,0.9), na.rm=T)
				lower_l_80[1,i] = as.numeric(quantiles[1])
				upper_l_80[1,i] = as.numeric(quantiles[2])
				traitMeanValue[1,i] = mean(traitValues[,i], na.rm=T)
				traitMedianValue[1,i] = median(traitValues[,i], na.rm=T)
			}
		tab = matrix(nrow=length(slicedTimes), ncol=6)
		colnames(tab) = c("time","median_value","95%HPD_lower_value","95%HPD_higher_value","80%HPD_lower_value","80%HPD_higher_value")	
		tab[,1] = slicedTimes; tab[,2] = traitMedianValue; tab[,3] = lower_l_95; tab[,4] = upper_l_95; tab[,5] = lower_l_80; tab[,6] = upper_l_80
		return(tab)
	}  

setwd(wd2)
registerDoMC(cores=1)
localTreesDirectory = "Minus_recombinant_regions_constraints_2_1000"
nberOfExtractionFiles = 1000
timeSlices = 100
slidingWindow = 10
traits = c(); hosts = c("Pan77","SC","GB")
measures = c("av_streak_width","chlorotic_area","deviation_from_rectangle","stunting","yellowness","whitness")
for (i in 1:length(measures))
	{
		for (j in 1:length(hosts)) traits = c(traits, paste(measures[i],hosts[j],sep="_"))
	}
for (i in 1:length(traits))
	{
		trait = traits[i]
		print(i); tab = traitEvolution(localTreesDirectory, nberOfExtractionFiles, timeSlices, slidingWindow, trait)
		write.csv(tab, paste0("Minus_recombinant_regions_constraints_values_2/Trait_values_",trait,".csv"), row.names=F, quote=F)
	}

setwd(wd2)
xLim = c(2007-224.3583, 2007) # MCC root height = 224.3583 with outgroup
xLim = c(2007-103.9772, 2007) # MCC root height = 103.9772 without outgroup
hosts = c("Pan77","SC","GB"); host_names = c("R","M","S"); cols = c("#1b9e77","#d95f02","#7570b3")
measures = c("av_streak_width","chlorotic_area","deviation_from_rectangle","stunting","yellowness","whitness")
measure_names = c("average streak width","chlorotic area","leaf deformation","leaf stunting","leaf yellowness","intensity of chlorosis")
pdf(paste("Minus_recombinant_regions_values_evolution_2.pdf",sep=""), width=6, height=5.5)
par(mfrow=c(2,2), mgp=c(1,0.35,0), oma=c(0.5,0.5,1.5,2), mar=c(3.3,3.3,0,0))
for (i in c(3,4,2,6))
	{
		tabs = list()
		for (j in 1:length(hosts))
			{
				trait = paste(measures[i],hosts[j],sep="_")
				tabs[[j]] = read.csv(paste0("Minus_recombinant_regions_constraints_values_2/Trait_values_",trait,".csv"), header=T)
				if (measures[i] == "stunting") tabs[[j]][,2:dim(tabs[[j]])[2]] = 1-tabs[[j]][,2:dim(tabs[[j]])[2]]
			}
		yMin = 10^10; yMax=-10^10
		for (j in 1:length(tabs))
			{
				if (yMin > min(tabs[[j]][,3])) yMin = min(tabs[[j]][,5])
				if (yMax < max(tabs[[j]][,4])) yMax = max(tabs[[j]][,6])
			}
		if (measures[i] == "av_streak_width") { yMin = 0.025; yMax = 0.06499137 }
		if (measures[i] == "stunting") { yMin = 0.25; yMax = 0.85 }
		yLim = c(yMin, yMax)
		for (j in 1:length(tabs))
			{
				if (j == 1)
					{
						plot(tabs[[j]][,1], tabs[[j]][,2], type="l", axes=F, ann=F, ylim=yLim, xlim=xLim, main=NULL, frame=F, col=NA)
						axis(side=1, lwd.tick=0.1, cex.axis=0.6, mgp=c(0,0.05,0), lwd=0, tck=-0.020, col.axis="gray30")
						axis(side=2, lwd.tick=0.1, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0, tck=-0.015, col.axis="gray30")
						title(xlab="time", cex.lab=0.7, mgp=c(1.0,0,0), col.lab="gray30")
						title(ylab=measure_names[i], cex.lab=0.7, mgp=c(1.3,0,0), col.lab="gray30")
						if (i == 6) legend(1950, 0.53, legend=rev(host_names), fill=rev(cols), border=0, text.col="gray30", bty="n", cex=0.7, y.intersp=1.2)
					}
				xx_l = c(tabs[[j]][,1], rev(tabs[[j]][,1]))
				yy_l = c(tabs[[j]][,5], rev(tabs[[j]][,6]))
				getOption("scipen"); opt = options("scipen"=20)
				polygon(xx_l, yy_l, col=rgb(col2rgb(cols[j])[1]/255,col2rgb(cols[j])[2]/255,col2rgb(cols[j])[3]/255,0.5), border=0)
				lines(tabs[[j]][,1], tabs[[j]][,2], lwd=1, col=cols[j])
				if (j == length(tabs)) box(lwd=0.1)
			}
	}
dev.off()

