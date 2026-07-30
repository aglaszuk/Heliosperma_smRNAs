# Title: "Detect differential targeting by small RNAs in Heliosperma samples grown in reciprocal transplantations (RT)"
# Author: Aglaia Szukala
# Short description: Script to perform differential targeting analyses (DT). Needs smRNAs counts table by region as input.

# Set project root (portable; works after cloning the repository)
.root_helper <- c("scripts/00_set_project_root.R", "../scripts/00_set_project_root.R")
source(.root_helper[file.exists(.root_helper)][1])
rm(.root_helper)

# Needs script 1_DE_functions.R 
source("scripts/1_DE_functions.R")

#Load Packages
library("RUVSeq")
library("edgeR")
library("HTSFilter")
library("RColorBrewer")
library("statmod")
library("gplots")
library("devtools")
library("ggpubr")
library("dplyr")

###################
# Load count data #
###################

x <- read.table(file = gzfile("data/smRNA_counts_RT.txt.gz"),
                header = TRUE,
                sep = "\t",
                stringsAsFactors = FALSE,
                check.names = FALSE
                )
head(x)
dim(x)
colnames(x)

####################
# Inspect raw data #
####################

# Make color vector by ecotype (shade) and ecotype pair (color)
specieColors <- c(rep("#A71E34",5), # alpine 1 down 
                  rep("#A71E34",5), # alpine 1 up 
                  rep("#00798C",5), # alpine 3 down
                  rep("#00798C",5), # alpine 3 up
                  rep("#EBADB5",5), # montane 1 down
                  rep("#EBADB5",5), # montane 1 up
                  rep("#B7E1D3",5), # montane 3 down
                  rep("#B7E1D3",2)  # montane 3 up
                  ) 

# Make symbols vector (ecotype)
geoNum <- c(
  rep(24, 5), # alpine 1 down: triangle
  rep(24, 5), # alpine 1 up
  rep(24, 5), # alpine 3 down
  rep(24, 5), # alpine 3 up
  rep(21, 5), # montane 1 down: circle
  rep(21, 5), # montane 1 up
  rep(21, 5), # montane 3 down
  rep(21, 2)  # montane 3 up
)

# Make filling vector (altitude) 
bgVec <- c(
  rep("#A71E34", 5), # alpine 1 down: filled
  rep(NA, 5),        # alpine 1 up: empty
  rep("#00798C", 5), # alpine 3 down: filled
  rep(NA, 5),        # alpine 3 up: empty
  rep("#EBADB5", 5), # montane 1 down: filled
  rep(NA, 5),        # montane 1 up: empty
  rep("#B7E1D3", 5), # montane 3 down: filled
  rep(NA, 2)         # montane 3 up: empty
)

# Plot PCA and save as pdf
pdf(file = "output/DT_RT/1_PCA/PCA_raw_counts.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
plotPCA(log2(cpm(x)+1),
        col=specieColors,
        pch = geoNum,
        bg = bgVec,
        k=2, 
        cex = 2, 
        cex.axis=1, 
        main = "PCA of untrimmed read counts", 
        labels = F
) 
dev.off()

# Check read counts distribution in raw data, as well as mean-variance relationship
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

# Define factors
snames <- colnames(rounded)
ecotype <- substr(snames, 1, nchar(snames) - 4) #Alpine vs Montane
ecotype
pair <- substr(snames, 2, nchar(snames) - 3) #ecotype pairs 1,3
pair
alt <- substr(snames, 4, nchar(snames) - 1) 
alt

# Create a new group variable
groupGLM <- interaction(ecotype, pair, alt)
groupGLM

# Trim data based on a min cpm treshhold in at least 3 individuals per group  
alpha <- 1
cpm_mat <- cpm(rounded)
min_samples <- 3

keep <- rowSums(
  sapply(levels(groupGLM), function(g) {
    rowSums(cpm_mat[, groupGLM == g, drop = FALSE] >= alpha) >= min_samples
  })
) > 0
y <- rounded[keep,]
dim(y)

# Plot data improvement upon trimming
par(mfrow=c(1,1))
pdf(file = "output/DT_RT/2_trimming/smRNA_counts_RT_trimming.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
hist(log(cpm(rounded)), 
     breaks= 400, 
     main=paste0("Distribution of log(cpm) in raw and trimmed data \nmean(cpm(count))>=",alpha),
     xlab= "Expression(log(cpm))", 
     col=rep("orange",400)
)
hist(log(cpm(y)), 
     breaks= 400, 
     col=rep("darkred",400),
     xlab= "Expression(log(cpm))", 
     add=T
)
dev.off()

##############################
# LIBRARY SIZE NORMALIZATION #
##############################

# Make a list-based data object, which contains a matrix of counts and a sample data-frame with samples and library size
y <- DGEList(counts=y,
             group=groupGLM
             ) 

# Look at variation in library size
par(mfrow=c(1,1))
pdf(file = "output/DT_RT/3_libSizeNormalization/LibrarySizeVariation.pdf",
    width = 8, # The width of the plot in inches
    height = 8)
with(y, plot(y$samples$group,log10(y$samples$lib.size),
             xlab="Groups", 
             ylab="Library size (log10)", 
             cex.axis=0.7,#x.axis=0.6, 
             main="Variation of the library size - Transplanted",
             col=(c("#A71E34","#EBADB5","#00798C","#B7E1D3")
                  )
             )
     ) 
dev.off()

# Compute normalization factors accounting for library size (seq depth)
y <- calcNormFactors(y, 
                     method ="TMM" # same used by Mimmi as ell
                     ) 
y$samples


# Save trimmed table of rounded counts
write.table(y$counts,file="output/DT_RT/2_trimming/smRNA_counts_RT_trimmed.txt",
            sep="\t",
            quote=F,
            col.names = T
            )

# Look at effect of normalization by sequencing depth
par(mfrow=c(2,2))
pdf(file = "output/DT_RT/3_libSizeNormalization/NormalizationByLibSize.pdf",   
    width = 8, # The width of the plot in inches
    height = 8)
makeRUVset(y$counts) # execute all graphs at once.  # execute all graphs at once. 
dev.off()

# Plot PCA of trimmed and normalized counts
par(mfrow=c(1,2))
pdf(file = "output/DT_RT/1_PCA/PCA_trimmed_normalized_counts.pdf",   
   width = 8, # The width of the plot in inches
   height = 8)
plotPCA(cpm(y), # Figure 1c in ms
        col=specieColors, 
        k=2, 
        cex = 1.5, 
        pch = geoNum, 
        bg = bgVec,
        cex.axis=1, 
        main = "PCA of normalized read counts", 
        labels = F
        ) 
dev.off()

##########################################################
# Detect DT regions using a Generalized log-Linear Model #
##########################################################

# Define a design matrix to describe the treatment conditions for the glm
MD <- model.matrix(~0+groupGLM)

# Define contrasts with Alpine ecotype as base level or native environment as base level
my.contrasts=makeContrasts(AvsM_1D=groupGLMM.1.D-groupGLMA.1.D, 
                           AvsM_3D=groupGLMM.3.D-groupGLMA.3.D,
                           AvsM_1U=groupGLMM.1.U-groupGLMA.1.U, 
                           AvsM_3U=groupGLMM.3.U-groupGLMA.3.U, 
                           MvsM1=groupGLMM.1.U-groupGLMM.1.D,
                           MvsM3=groupGLMM.3.U-groupGLMM.3.D,
                           AvsA1=groupGLMA.1.D-groupGLMA.1.U,
                           AvsA3=groupGLMA.3.D-groupGLMA.3.U,
                           levels=MD)
my.contrasts

# Use the EdgeR method to detect DT regions
outdir <- "output/DT_RT/4_DT_edgeR/"

# Estimate pairwise dispersion
DGE=estimateDisp(y, 
                 design = MD, 
                 robust = T
                 )

# Dispersion (sd) around the mean expression level
sqrt(DGE$common.dispersion)

# Fit a quasi-likelihood negative binomial generalized log-linear model to count data
GLMql = glmQLFit(DGE,
                 design = MD,
                 robust=T
                 )

# Perform region-wise statistical tests for a given contrast and identify DT regions 
# FDR threshold for significance: 0.05
par(mfrow=c(2,2))
findDTregions(GLMql,my.contrasts[,1], "AvsM_1D")
findDTregions(GLMql,my.contrasts[,2], "AvsM_3D")
findDTregions(GLMql,my.contrasts[,3], "AvsM_1U")
findDTregions(GLMql,my.contrasts[,4], "AvsM_3U")
findDTregions(GLMql,my.contrasts[,5], "MvsM1")
findDTregions(GLMql,my.contrasts[,6], "MvsM3")
findDTregions(GLMql,my.contrasts[,7], "AvsA1")
findDTregions(GLMql,my.contrasts[,8], "AvsA3")