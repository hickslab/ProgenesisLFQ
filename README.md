# *`ProgenesisLFQ`*

**This repository contains functions to process Progenesis QI for proteomics v2.0 data. Download a workflow file depending on your application and open it in RStudio to edit/run. Find '###' in the workflows and update with your particular filename/path before running. The following data files are needed:**

* Protein measurements:
  + [Global proteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Global_Workflow.R)
  
* Peptide measurements, protein measurements, and FASTA database:
  + [Phosphoproteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Phospho_Workflow.R)
  + [Redox proteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Redox_Workflow.R)


# *`StatLFQ`*

**Statistical workflow for differential analysis. Takes ProgenesisLFQ output in simple format (variable column followed by abundance values). The following tests are available:**

* [Pairwise *t*-test](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/StatLFQ.R)
  
* [One-way ANOVA](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/StatLFQ.R)


# *`AnnotateLFQ`*

* [Remove PACid](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)

* [KEGG](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)

* [GO](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)


# *`PlotLFQ`*

* [PCA](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/PlotLFQ.R)