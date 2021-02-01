gsea.compiler<-function(posFile,negFile, sep = ",", outputname = "compiledGSEA.csv") {
  pos <- read.table(posFile, sep = sep, header = TRUE)
  neg <- read.table(negFile, sep = sep, header = TRUE)
  compiled <- rbind(pos,neg)
  write.csv(compiled, outputname, row.names = F)
}
