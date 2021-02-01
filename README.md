# GSEA-Plotting-Tools
Functions to visualize results from GSEA. 

GSEA scatter is used to create a GSEA plot (example in repository) with one or two sample files. If plotting two sample GSEA data sets, both analyses must contain the same gene sets. The genesetlist file (example in repository) contains the names of the genesets in the GSEA analysis files that you want plotted. The sample filenames will be used as the sample labels in the final plot.

```
usage:
genesetlist = file containing the list of gene sets of interests
gseacsvs =  vector containing the string/strings of gsea data files to be plotted
graphname = output pdf graph name
width = width of output graph
height = height of output graph

gsea.scatter("c2_enriched_GCN2_ER_AA_4pathways.txt", c("6H_Gln_starved.csv","6H_Tunicamycin.csv"), graphname = "sampleGSEAplot.pdf")
```


GSEA compiler is a basic function to combine the upregulated and downregulated GSEA results files into a single .csv file that can be used with the GSEA scatter function.
```
usage:
posFile = upregulated GSEA results file
negFile = downregulated GSEA results file
sep = delimiter for the input results files (default ",")
outputname = output filename

gsea.compiler(posFile = "posfile.csv", negFile = "negfile.csv", sep = ",", outputname = "compiledFile.csv")
```
