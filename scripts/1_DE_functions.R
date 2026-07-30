#Compare methods to trimm read counts using cpm
trimMetComp <- function(dat,groups,minCount){ 
  mm <- model.matrix(~0+groups)
  #method 1
  k <- filterByExpr(dat, design = mm, min.count=minCount)
  out1 <- dat[k,]
  print(paste0("Number of genes before trimming:",nrow(dat)))
  print(paste0("Number of genes after trimming using build in function:",nrow(out1)))
  out1 <- DGEList(counts=out1,group=groups)
  out1 <- calcNormFactors(out1, method ="TMM")
  
  #method 2
  keep <- rowSums(cpm(dat)>=minCount) >= 4 #number of biological replicates
  out2 <- x[keep,]
  print(paste0("Number of genes after trimming using min count in min 4 samples:",nrow(out2)))
  out2 <- DGEList(counts=out2,group=groups)
  out2 <- calcNormFactors(out2, method ="TMM")
  
  #method 3
  out3=dat[apply(cpm(dat),1,function(x){!(mean(x)<minCount)}),]
  print(paste0("Number of genes after trimming using min mean:",nrow(out3)))
  out3 <- DGEList(counts=out3,group=groupGLM)
  out3 <- calcNormFactors(out3, method ="TMM")
  
  #plot 
  par(mfrow=c(1,3))
  voom(out1, mm, plot = T)
  voom(out2, mm, plot = T)
  voom(out3, mm, plot = T)
  
  #plot raw and trimmed count distributions
  par(mfrow=c(1,1))
  hist(log(cpm(dat)), breaks= 1000, col=rep("antiquewhite2",1000),
       main="Distribution of log(cpm) in raw data",
       xlab= "Expression(log(cpm))")
  hist(log(cpm(out1)), breaks= 1000,
       main=paste0("Distribution of log(cpm) after trimming (mean(cpm)>",alpha,")"),
       xlab= "Expression(log(cpm))",col=rep("darkgoldenrod2",1000),
       xlim=c(-7,10), add=T)
  hist(log(cpm(out2)), breaks= 1000,
       main=paste0("Distribution of log(cpm) after trimming (mean(cpm)>",alpha,")"),
       xlab= "Expression(log(cpm))",col=rep("coral2",1000),
       xlim=c(-7,10), add=T)
  hist(log(cpm(out3)), breaks= 1000,
       main=paste0("Distribution of log(cpm) after trimming (mean(cpm)>",alpha,")"),
       xlab= "Expression(log(cpm))",col=rep("darkcyan",1000),
       xlim=c(-7,10), add=T)
}

#Plot data before and after removing unwanted variance (normalization step)
makeRUVset <- function (dat){
  plotRLE(as.matrix(dat),
          outline=FALSE,
          ylim=c(-4,4),cex.axis=0.5, 
          col=specieColors,
          main = "Unwanted variance not removed")
  plotPCA(as.matrix(dat), 
          col=specieColors, 
          k=2, cex=1, pch = geoNum, labels = F,
          cex.axis=1, main= "PCA - Unwanted variance not removed")
  # RLE-plot -> Visualizing unwanted variation in high dimensional data
  uq = betweenLaneNormalization(as.matrix(dat), which = "full")
  set = newSeqExpressionSet(uq)
  plotRLE(set,
          outline=FALSE,
          ylim=c(-4,4),col=specieColors, 
          cex.axis=0.5, main = "Unwanted variance removed")
  
  plotPCA(set,
          col=specieColors, k=2, cex = 1, 
          pch = geoNum, cex.axis=1, 
          main = "PCA - Unwanted variance removed", labels = F) 
  return(set)
}

#Find DT regions using voom function and limma package
findDTregVoom <- function(dat, model, contr, name){
  #apply weights to samples such that outlier samples are downweighted
  v <- voom(dat, model, normalize.method="none", plot = F) #voomWithQualityWeights cam be used alterntively
  #lmFit fits a linear model using weighted least squares for each gene
  fit <- lmFit(v, model) 
  tmp <- contrasts.fit(fit, contr)
  #robust=TRUE setting to leverage the quality weights such that the analysis is robust to outliers
  tmp <- eBayes(tmp, robust=TRUE) 
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf, adjust.method="BH", p.value=0.05)
  nup <- length(rownames(top.table[top.table$logFC>0 & top.table$adj.P.Val<0.05,]))
  ndw <- length(rownames(top.table[top.table$logFC<0 & top.table$adj.P.Val<0.05,]))
  print(paste0("Number of uptargeted regions: ",nup))
  print(paste0("Number of downtargeted regions: ",ndw))
  
  write.table(top.table[top.table$adj.P.Val<0.05,], 
              file = paste0(outdir,name,"_DTreg.txt"), 
              sep = "\t", quote = F)
}

#Find DT regions using voom with quality weights function and limma package
findDTregVoomWQW <- function(dat, model, contr, name){
  #apply weights to samples such that outlier samples are downweighted
  v <- voomWithQualityWeights(dat, model, normalize.method="none", plot = F) #voomWithQualityWeights cam be used alterntively
  #lmFit fits a linear model using weighted least squares for each gene
  fit <- lmFit(v, model) 
  tmp <- contrasts.fit(fit, contr)
  #robust=TRUE setting to leverage the quality weights such that the analysis is robust to outliers
  tmp <- eBayes(tmp, robust=TRUE) 
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf, adjust.method="BH", p.value=0.05)
  nup <- length(rownames(top.table[top.table$logFC>0 & top.table$adj.P.Val<0.05,]))
  ndw <- length(rownames(top.table[top.table$logFC<0 & top.table$adj.P.Val<0.05,]))
  print(paste0("Number of uptargeted regions: ",nup))
  print(paste0("Number of downtargeted regions: ",ndw))
  
  write.table(top.table[top.table$adj.P.Val<0.05,], 
              file = paste0(outdir,name,"_DTreg.txt"), 
              sep = "\t", quote = F)
  
  #volcano plot
  top.table$de <- "no"
  top.table$de[top.table$logFC>2] <- "up"
  top.table$de[top.table$logFC<(-2)] <- "dw"
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("dw", "up", "no")
  # top.table$delabel <- NA
  # top.table$delabel[top.table$de != "no"] <- rownames(top.table)[top.table$de != "no"]
  ggplot(data=top.table, aes(x=logFC, y=-log10(adj.P.Val), col=de)) +
    geom_point() +
    theme_minimal() +
    geom_vline(xintercept=c(-2, 2), col="red") +
    geom_hline(yintercept=-log10(0.05), col="red") +
    scale_colour_manual(values = mycolors)
  #geom_text()
  #ma plot
  ggplot(data=top.table, aes(x=AveExpr, y=logFC, col=de)) + #AveExpr is log2 transformed!
    geom_point() +
    theme_minimal() +
    xlab("log2 Average expression") +
    geom_hline(yintercept=c(-2, 2), col="red") +
    scale_colour_manual(values = mycolors)
}

#Find DT regions using edgeR and the Quasi-Likelihood framework
findDTregions <- function(glm, contrast, name){
  qlt <- glmQLFTest(glm,contrast = contrast)
  de.genes = topTags(qlt,p.value = 0.05, adjust.method = "BH",n=dim(y)[1])
  write.table(de.genes$table[de.genes$table$logFC > 0,], file=paste0(outdir,name,"_DTRs_qlTest_poslogFC.txt"),
              quote = F, sep = "\t")
  write.table(de.genes$table[de.genes$table$logFC < 0,], file=paste0(outdir,name,"_DTRs_qlTest_neglogFC.txt"),
              quote = F, sep = "\t")
  
  #Number of up-targeted regions in comparison
  nUp <- length(rownames(de.genes$table[de.genes$table$FDR <= 0.05 & de.genes$table$logFC>0,]))
  print(paste0("Number of uptargeted regions: ",nUp))
  #Number of down-targeted regions in comparison
  nDw <- length(rownames(de.genes$table[de.genes$table$FDR <= 0.05 & de.genes$table$logFC<(-0),]))
  print(paste0("Number of downtargeted regions: ",nDw))
  
  top.table <- de.genes$table
  top.table$de <- "no"
  top.table$de[top.table$FDR < 0.05] <- "yes"
  mycolors <- c("firebrick3","grey")#"blue", "red",
  names(mycolors) <- c("yes","no")
  ggplot(data=top.table, aes(x=logCPM, y=logFC, col=de)) + #AveExpr is log2 transformed!
    geom_point() +
    theme_minimal() +
    xlab("log CPM") +
    geom_hline(yintercept=c(-0.6, 0.6), col="red") +
    scale_colour_manual(values = mycolors)
  plotSmear(qlt,de.tags = rownames(de.genes$table)[which(de.genes$table$FDR<0.05)],
            xlab="Average logCPM ", ylab="logFC",
            main=paste("Differential Targeted Regions in ",name,"\nUptargeted:", nUp, ", Downtargeted:", nDw),
            cex.sub=0.7,
            xlim=c(-2,16),
            ylim=c(-11,13))
  points(qlt$table$logCPM,
         qlt$table$logFC,
         col="gray18",pch=16, cex=0.7)
  points(qlt$table$logCPM[which(rownames(qlt$table) %in% rownames(de.genes$table))],
         qlt$table$logFC[which(rownames(qlt$table) %in% rownames(de.genes$table))],
         col="black", bg="firebrick3", pch=21, cex=0.7)
}

#Find DT regions using edgeR and the Likelihood-Ratio test
findDTgenesLRT <- function(glm,contrast, name){
  lrt <- glmLRT(glm,contrast = contrast)
  de.genes = topTags(lrt,p.value = 0.05, adjust.method = "BH",n=dim(y)[1])
  write.table(de.genes$table, file=paste0(outdir,name,"_DEregions_lrt.txt"),
              quote = F, sep = "\t")
  
  #Number of up-targeted regions in comparison
  nUp <- length(rownames(de.genes$table[de.genes$table$FDR <= 0.05 & de.genes$table$logFC>0,]))
  print(paste0("Number of uptargeted regions: ",nUp))
  #Number of down-targeted regions in comparison
  nDw <- length(rownames(de.genes$table[de.genes$table$FDR <= 0.05 & de.genes$table$logFC<0,]))
  print(paste0("Number of downtargeted regions: ",nDw))
  
  plotSmear(lrt,de.tags = rownames(de.genes$table)[which(de.genes$table$FDR<0.05)],
            xlab="Average logCPM ", ylab="logFC", 
            main=paste("Differential Targeted Regions in ",name,"\nUptargeted:", nUp, ", Downtargeted:", nDw), 
            cex.sub=0.7,
            xlim=c(-2,16),
            ylim=c(-11,11))
}

#Find DT regions dected by both voom and edgeR and plot intersection
DTR_intersect <- function(y,GLM, contr, color, name){
  #find DEGs using voom
  v <- voom(y, MD1, normalize.method="none", plot = F)
  fit <- lmFit(v, MD1)
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp, robust=TRUE)
  top.table <- topTable(tmp, sort.by="P", adjust.method="BH", p.value=0.05, n=Inf)
  #find DEGs using edgeR
  lrt <- glmQLFTest(GLM,contrast=contr)
  de.genes = topTags(lrt,p.value = 0.05, adjust.method = "BH",n=dim(y)[1])
  res <- de.genes$table
  #find intersection
  res_intersect <- res[rownames(res) %in% rownames(top.table),]
  write.table(res_intersect, file = paste0(outdir,name,".txt"),
              quote = F, sep = "\t")
  
  print(paste0("The number of DTR identified by voom: ", dim(top.table)[1]))
  print(paste0("The number of DTR identified by edgeR: ", dim(res)[1]))
  print(paste0("The number of DTR identified by both methods is: ", dim(res_intersect)[1]))
  
  nUp <- length(rownames(res_intersect[res_intersect$FDR <= 0.05 & res_intersect$logFC>0,]))
  nDw <- length(rownames(res_intersect[res_intersect$FDR <= 0.05 & res_intersect$logFC<0,]))
  
  #volcano plot
  top.table$de <- "no_intersect"
  top.table$de[which(rownames(top.table) %in% rownames(res))] <- "shared"
  # top.table$de[top.table$logFC>2] <- "up"
  # top.table$de[top.table$logFC<(-2)] <- "dw"
  mycolors <- c("firebrick3","grey")#"blue", "red", 
  names(mycolors) <- c("shared","no_intersect")#"dw", "up", 
  #ma plot
  ggplot(data=top.table, aes(x=AveExpr, y=logFC, col=de)) + #AveExpr is log2 transformed!
    geom_point() +
    theme_minimal() +
    xlab("log2 Average expression") +
    # ggtitle(paste0("The number of DTR identified by voom: ", dim(top.table)[1],
    #      "\nThe number of DTR identified by edgeR: ", dim(res)[1]),
    #      "\nThe number of DTR identified by both methods is: ", dim(res_intersect)[1])+
    geom_hline(yintercept=c(-0.6, 0.6), col="red") +
    scale_colour_manual(values = mycolors)
  # plotSmear(lrt,
  #           de.tags = rownames(lrt$table)[which(rownames(lrt$table) %in% rownames(res_intersect))],
  #           xlab="Average logCPM ", ylab="logFC",
  #           main=paste(name,"\nUptargeted:", nUp, ", Downtargeted:", nDw),
  #           cex.sub=0.7,
  #           xlim=c(0,16),
  #           ylim=c(-10,10))
  # points(lrt$table$logCPM,
  #        lrt$table$logFC,
  #        col="antiquewhite3",pch=16)
  # points(lrt$table$logCPM[which(rownames(lrt$table) %in% rownames(res))],
  #        lrt$table$logFC[which(rownames(lrt$table) %in% rownames(res))],
  #        col="antiquewhite4", pch=16)
  # points(lrt$table$logCPM[which(rownames(lrt$table) %in% rownames(res_intersect))],
  #        lrt$table$logFC[which(rownames(lrt$table) %in% rownames(res_intersect))],
  #        col="black",pch=21, bg=color)
}