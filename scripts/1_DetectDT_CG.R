# Title: "Detect differential targeting (DT) by small RNAs in a common garden (CG) experimental setup"
# Author: Aglaia Szukala
# Short description: Script to perform differential targeting analyses (DT). Needs a table of smRNAs counts by region as input.
# Needs script 1_DE_functions.R 

# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Load functions script
source("scripts/1_DE_functions.R")

#Load Libraries
library("RUVSeq") # Remove Unwanted Variation from RNA-Seq Data
library("edgeR") # Perform DT analysis
library("HTSFilter")
library("RColorBrewer")
library("statmod")
library("gplots")
library("devtools")
library("ggpubr")
#library("VennDiagram") 
#library("pheatmap")

###################
# Load count data #
###################

x <- read.table(file = gzfile("data/smRNA_counts_CG.txt.gz"))
head(x)

# Define factors
snames <- colnames(x)
ecotype <- substr(snames, 1, nchar(snames) - 2) #Alpine vs Montane
pair <- substr(snames, 2, nchar(snames) - 1) #ecotype pairs 1,3,4 and 5

# Define colors and symbols for samples mirroring localities (colors) and ecotypes (shades, symbols)
specieColors <- c(rep("#A71E34",3), 
                  rep("#00798C",2),
                  rep("#F59700",3),
                  rep("#6E5DAC",3),
                  rep("#EBADB5",4), 
                  rep("#B7E1D3",3),
                  rep("#FFD085",3),
                  rep("#CEACEC",2)
                  )

geoNum <- c(rep(17,11),
            rep(19,12)
)

########################
#  RAW DATA INSPECTION #
########################

# Inspect raw data by means of PCA and save plot as pdf
pdf(file = "output/DT_CG/1_PCA/PCA_raw_counts.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
plotPCA(log2(cpm(x)+1),
        col=specieColors, k=2, cex = 1, 
        #pch = geoNum, 
        cex.axis=1, 
        labels = T
        ) 
dev.off()

# Inspect read counts distribution in raw data
par(mfrow=c(1,2))
hist(log(cpm(x)), 
     breaks= 1000, 
     main="Distribution of log(cpm) in raw data",
     xlab= "Expression(log(cpm))"
     )

plot(apply(log(cpm(x)), 1, mean), 
     apply(log(cpm(x)), 1, var),
     xlab="Mean expression (log(cpm))", 
     ylab="Variance of expression log(cpm)",
     main="Raw - Expression variance"
     )

##################
#  DATA TRIMMING #
##################

# Round count values to obtain integers
rounded <- round(x, 
                 digits = 0
                 )

# Trim regions using a min threshold of ca 5-10 reads in at least 3 samples
alpha <- 1
keep <- rowSums(cpm(rounded)>=alpha) >= 3 
y <- rounded[keep,]

# Save trimmed table of rounded counts
write.table(y,file="output/DT_CG/2_trimming/smRNA_counts_CG_trimmed.txt",
            sep="\t",
            quote=F,
            col.names = T
)

# Check amount of retained regions after trimming
dim(x) # 731619 before trimming
dim(y) # 286427 after trimming

# Plot distribution of CPM before and after trimming and save as pdf
par(mfrow=c(1,1))
pdf(file = "output/DT_CG/2_trimming/smRNA_counts_CG_trimming.pdf",   
    width = 8, 
    height = 8)
hist(log(cpm(rounded)), 
     breaks= 400, 
     main=paste0("Distribution of log(cpm) in raw (blue) and trimmed (orange) data \nmean(cpm(count))>=",alpha),
     xlab= "Expression(log(cpm))", 
     col=rep("blue",400)
     )
hist(log(cpm(y)), 
     breaks= 400, 
     col=rep("orange",400),
     xlab= "Expression(log(cpm))", 
     add=T
     )
dev.off()

# Plot PCA of trimmed counts
par(mfrow=c(1,2))
pdf(file = "output/DT_CG/1_PCA/PCA_trimmed_normalized_counts.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
plotPCA(cpm(y), # Figure 1c in ms
        col=specieColors, 
        k=2, 
        cex = 2, 
        pch = geoNum, 
        cex.axis=1, 
        main = "PCA of normalized read counts", 
        labels = F
) 
dev.off()

##############################
# LIBRARY SIZE NORMALIZATION #
##############################

# Create a new group variable
groupGLM <- interaction(ecotype, pair)
groupGLM

# Make a list-based data object, which contains a matrix of counts and a sample data-frame with samples and library size
y <- DGEList(counts=y,
             group=groupGLM
             ) 

# Look at variation in library size
par(mfrow=c(1,1))
pdf(file = "output/DT_CG/3_libSizeNormalization/LibrarySizeVariation.pdf",
    width = 8, # The width of the plot in inches
    height = 8)
with(y, plot(y$samples$group,log10(y$samples$lib.size),
             xlab="Populations", 
             ylab="Library size (log10)", 
             cex.axis=0.7,#x.axis=0.6, 
             main="Variation of the library size - Transplanted",
             col=(c("#A71E34","#EBADB5",
                    "#00798C","#B7E1D3",
                    "#F59700","#FFD085",
                    "#6E5DAC","#CEACEC")
                  )
             )
     ) 
dev.off()

# Look at the effect of normalization by library size
par(mfrow=c(2,2))
pdf(file = "output/DT_CG/3_libSizeNormalization/NormalizationByLibSize.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
makeRUVset(y$counts)
dev.off()

# Compute normalization factors accounting for library size
y <- calcNormFactors(y, 
                     method ="TMM"
                     ) 

# Inspect the new object
y$samples # normalization factors

##########################################################
# Detect DT regions using a Generalized log-Linear Model #
##########################################################

# Define a design matrix to describe the treatment conditions for the glm removing the intercept
MD <- model.matrix(~0+groupGLM)

# Define contrasts with Alpine ecotype as base level
my.contrasts=makeContrasts(AvsM_1=groupGLMM.1-groupGLMA.1, 
                           AvsM_3=groupGLMM.3-groupGLMA.3,
                           AvsM_4=groupGLMM.4-groupGLMA.4, 
                           AvsM_5=groupGLMM.5-groupGLMA.5, 
                           levels=MD
                           )
my.contrasts

# Use the EdgeR method to detect DT regions
# Define output directory
outdir <- "output/DT_CG/4_DT_edgeR/"

# First estimate pairwise dispersion based on model design
DGE=estimateDisp(y, 
                 design = MD, 
                 robust = T
                 )

# Dispersion (sd) around the mean expression level

# Fit a quasi-likelihood negative binomial generalized log-linear model to count data
GLMql = glmQLFit(DGE,
                 design = MD,
                 robust=T
                 )

# Perform region-wise statistical tests for a given contrast and identify DT regions 
# FDR threshold for significance: 0.05
par(mfrow=c(2,2))
findDTregions(GLMql,my.contrasts[,1], "MvsA1")
findDTregions(GLMql,my.contrasts[,2], "MvsA3")
findDTregions(GLMql,my.contrasts[,3], "MvsA4")
findDTregions(GLMql,my.contrasts[,4], "MvsA5")

