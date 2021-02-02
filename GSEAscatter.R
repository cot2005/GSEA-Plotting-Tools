library(ggplot2)

gsea.scatter<-function(genesetlist, gseacsvs, graphname = "GSEAplot.pdf", width = 8, height = 5) {
  targetgeneset <- readLines(genesetlist)
  samplenames <- gsub(".csv", "", gseacsvs)
  numsamples <- length(gseacsvs)
  for (i in 1:numsamples) {
    tempgsea <- read.csv(gseacsvs[i], stringsAsFactors = F)
    tempgsea <- tempgsea[match(targetgeneset, tempgsea$NAME),]   #matches genesets
    tempgsea$Samplename <- samplenames[i]
    #takes the samplename, geneset name, NES, FDR, leading edge
    tempgsea <- subset(tempgsea, select=c(12,1,6,8,11))
    if (i == 1) {
      gseadf <- tempgsea
    } else {
      gseadf <- rbind(gseadf, tempgsea)
    }
  }
  gseadf$LEADING.EDGE <- substr(gseadf$LEADING.EDGE,6,8)   #grabs tag percentage
  gseadf$LEADING.EDGE <- gsub("%", "", gseadf$LEADING.EDGE)
  gseadf$LEADING.EDGE <- gsub(",", "", gseadf$LEADING.EDGE)
  gseadf$LEADING.EDGE <- as.numeric(gseadf$LEADING.EDGE)/100
  
  colnames(gseadf) <- c("Sample","GeneSet", "NES","FDR","GeneSetProportion")
  gseadf$GeneSet <- factor(gseadf$GeneSet, levels = gseadf$GeneSet[length(targetgeneset):1])
  narows <- which(is.na(gseadf$GeneSet == TRUE))   #removes NAs from table
  if (length(narows) > 0) {
    gseadf <- gseadf[-narows,]
  }
  
  # calculates graphing parameters
  nesmax <- max(gseadf$NES) * 1.05
  if (nesmax < 0) {
    nesmax <- max(gseadf$NES) * 0.95
  } 
  nesmin <- min(gseadf$NES) * 0.95
  if (nesmin < 0 ) {
    nesmin <- min(gseadf$NES) * 1.05
  }
  if (numsamples == 1) {
    numshapes <- 5
  } else {
    numshapes <- 3
  }
  propmax <- floor(max(gseadf$GeneSetProportion * 10)) + 0.5
  if (propmax > max(gseadf$GeneSetProportion * 10)) {
    propmax <- floor(max(gseadf$GeneSetProportion * 10))
  }
  propmin <- ceiling(min(gseadf$GeneSetProportion * 10)) - 0.5
  if (propmin < min(gseadf$GeneSetProportion * 10)) {
    propmin <- ceiling(min(gseadf$GeneSetProportion * 10))
  }
  proprange <- signif(seq(propmin, propmax, (propmax - propmin)/(numshapes - 1)), digits = 2)
  proprange <- proprange/10
  shapelist <- rep(21, numshapes)
  shapekey <- 21
  if (numsamples == 2) {
    shapelist <- append(shapelist, rep(24, numshapes))
    shapekey <- append(shapekey, 24)
  }
  
  #plots data. NES is x axis, gene set in alphabetical on Y axis, colored is FDR, Shape is sample, size is geneset proportion
  g <- ggplot(gseadf, aes(x = NES, y = GeneSet, shape = Sample, size = GeneSetProportion, color = FDR))
  g <- g + geom_point() + scale_color_gradientn(colors = c("red", "blue2", "azure3"), breaks = c(0.05,0.25,0.5,0.75,1) ,values = scales::rescale(c(0,.05,.5)), guide = guide_colorbar(reverse = TRUE), limits = c(0,1)) + 
    scale_shape_manual(values = c(16,17)) + theme_bw() + coord_cartesian(xlim = c(nesmin, nesmax)) +
    guides(size = guide_legend(override.aes = list(shape = shapelist)), shape = guide_legend(override.aes = list(shape = shapekey, size = 5))) + 
    labs(title = "GSEA Comparison", x = "Normalized Enrichment Score (NES)") + scale_size_continuous(name = "Enriched Geneset\nProportion", breaks = rep(proprange,numsamples), range=c(2,8))
  print(g)
  ggsave(graphname, width = width, height = height)
}
