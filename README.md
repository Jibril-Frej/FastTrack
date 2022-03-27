# FastTrack

Tutorial to perform Differential Gene Expression Analysis with R and RStudio.

## Install (Windows 10 only)
Download and install R : [https://cran.r-project.org/bin/windows/base/](https://cran.r-project.org/bin/windows/base/)

Download and install RStudio : [https://www.rstudio.com/products/rstudio/download/#download](https://www.rstudio.com/products/rstudio/download/#download)

In RStudio Console type the following commands to install the libraries needed by FastTrack:

Install Bioconductor 
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.13")
```

Install DESeq2 and hgu133plus2.db from Bioconductor
```R
BiocManager::install("DESeq2")
BiocManager::install("hgu133plus2.db")
```

Install all other libraries
```R
install.packages("pheatmap")
install.packages("ggplot2")
install.packages("readxl")
install.packages("writexl")
install.packages("tibble")
install.packages("scales")
install.packages("viridis")
install.packages("tidyr")
```

## Usage
Download the repository and open the file UserGuide.nb.html (or AdvancedUserGuide.nb.html) with your browser and follow the instructions.


## Sources
Part of the code comes from : 

[https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/](https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/)

[https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html](https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html)
