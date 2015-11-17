# Enrich_shiny
R/shiny app for enrichment analysis

This app requires the following packages : shiny, org.Hs.eg.db, org.Mm.eg.db, GO.db, allez, EACI

To install shiny and db packages, in R run:

install.packages("shiny")

source("http://bioconductor.org/biocLite.R")

biocLite("org.Hs.eg.db”)

biocLite("org.Mm.eg.db”)

biocLite("GO.db")

allez and EACI could be found at https://github.com/lengning/Enrichment-test/tree/master/pkgs

To install these two packages, in bash run 

R CMD INSTALL allez.tar.gz

R CMD INSTALL EACI_0.0.1.tar.gz


The input file should contain gene scores. The score could be either continuous or binary (1 for on and 0 for off).
Note for EASE method, only binary scores could be used as input.
The input file should contain only one column. Row names should be gene names. Each entry represents loading (or score) 
of that particular gene.
Currently the program takes csv files or tab delimited file.
The input file will be treated as a tab delimited file if the suffix is not '.csv'.
Example input files could be found at https://github.com/lengning/Enrich_shiny/tree/master/example_input   (Exampleloadings)

The second input button could be used to input the marker lists of interest. It could be csv or tab delimited file. Each row represents a marker list. 
Row names are the list names. The number of columns should be N, in which N is the number of genes in the longest list. 
For lists with length M - N, the M+1, ..., N's column in that row should be filled with "" or " ". If the local marker list
is not specified, only GO terms will be considered. 
Example input files could be found at https://github.com/lengning/Enrich_shiny/tree/master/example_input   (MarkerLists)

The threshold to filter out small(large) sets defaults to 5 (500). Sets that are smaller (larger) than the threshold are not considered.

The annotation could be either 'human' or 'mouse'. The p value cutoff defaults to 0.1.

The user could also choose whether one-tailed p value should be calculated. If the input values are absolute values, then it should be specified as "one-tailed". Otherwise it should be specified as "two-tailed", in which case two-tailed p values will be calculated. This option is only valid in the allez implementation (not in EACI or EASE).

The user could use the 'output folder' button to choose the output directory. If it is not selected, home directory (~/) will be used
as the output directory. 

Outputs:

XX_allsets.txt: enrichment results of all gene sets (GO sets + local sets). 

XX_allsets_significant.txt: significant terms; The summary statistics of these sets are the same as those in the XX_allsets.txt

XX_localsets.txt: enrichment results of only local sets. The summary statistics of these sets are the same as those in the XX_allsets.txt

XX_localsets_significant.txt: significant terms in local lists


The output files contain GO term (NA for local sets), set p value, adjusted p value, z score, set size, set mean and set sd. Sets are sorted by p value. Sets with large absoloute z scores are expected to have small p values.

