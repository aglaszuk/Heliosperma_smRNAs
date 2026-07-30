topGoAnalysis <- function(ontology, allGenes, gene2GO, nbNodes){
  GOdata <- new("topGOdata", ontology = ontology, allGenes = allGenes, annot = annFUN.gene2GO, gene2GO = gene2GO)
  resultWeight01 <- runTest(GOdata, statistic = "fisher")
  allRes <- GenTable(GOdata, weight01_pval=resultWeight01, orderBy = "weight01", ranksOf = "weight01",topNodes = nbNodes)
  allRes <- cbind(allRes,"BP")
  colnames(allRes) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  return(allRes)
}

change_names <- function(data, name_list){
  colnames(data) <- name_list
  return(data)
}

rename <- function(table, geneNames){
  names(table) <- geneNames
  return(table)
}

attach_enriched_go_genes <- function(enriched_go_with_my_genes){
  enriched_go_with_my_genes.list = c()
  for (i in 1:length(enriched_go_with_my_genes)){
    enriched_go_with_my_genes.list = c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
  }
  return(enriched_go_with_my_genes.list)
}
circle_dat <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- logFC[s:e]
    #value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    zsc <- c(zsc, sum(value, na.rm = F) / sqrt(count[c])) #### HERE : na.rm = TRUE, takes all the genes into account !
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

plot_GO <- function(data_GO){
  #read in data matrix
  data <- read.table(data_GO, sep = "\t")
  colnames(data) <- c("pair","ecotype","category","ID","term","count","adj_pval","zscore")
  # define colors of bars by zscores
  ggplot(data, aes(x=reorder(term,zscore),y=-log10(adj_pval)))+ #x=reorder(term,zscore),y=-log10(adj_pval)
    geom_col(aes(fill = zscore)) +
    scale_fill_gradient(low="blueviolet",high="goldenrod1", 
                        labels=c(paste0("Underexpressed","\n","in L"),
                                 paste0("Overexpressed","\n", "in L")),
                        breaks=c(min(data$zscore),max(data$zscore)),
                        limits=c(min(data$zscore),max(data$zscore))) +
    geom_text(aes(label=count), 
              hjust=-1, size=5) +
    xlab("")+
    theme_bw()+
    coord_flip(ylim = c(0, 5.3)) +
    #scale_y_reverse() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    #coord_cartesian(ylim = c(0, 4.6))+
    theme(axis.text.x = element_text(
      angle = 0,
      hjust = 0.3,
      vjust = 0.5))+
    theme(text=element_text(family="Arial"))+
    theme(axis.text=element_text(size=15), axis.text.y = element_text(size=15)) +
    #       axis.title=element_text(size=14,face="bold"))
    theme(legend.text=element_text(size=13),
          legend.title = element_text(size=15, vjust=2.1),
          axis.title.x = element_text(size=15))
}