# Puig et al., 2023: bulk and single-cell RNAseq analysis scripts
 
This repository simply captures the analysis and the rationale behind the approaches taken.

## Description

Scripts #1-4 deal with the bulk RNAseq analysis of data generated by Jerome Korzelius and colleagues. Analysis pipeline written by Joaquín de Navascués, based on earlier work by Aleix Puig-Barbé.

Scripts #5-7 deal with the analysis of scRNAseq data already published by the Perrimon lab and the Fly Cell Atlas consortium. The pipeline was written by Vinicius Dias Nirello, and adapted for sharing by Joaquín de Navascués.

All scripts are designed so that the repository can be cloned and run straight away. The `librarian` package should take care of installing whatever packages are needed, and all resources are either available or downloaded programmatically. The tests have only been done in RStudio, though.

## In-progress bits

- There is a mismatch between the Rmd and HTML versions of script #1 - the HTML relies on a previous version, where the read count data per gene per sample was stored as a zip file in the repo. Now the Rmd file relies on reading the data directly from GEO. This will be tested and knitted into HTML as soon as the GEO data is released (should be by 12/06/2023).
- Scripts #5-7 are not ready yet (07/06/2023), but will be made available soon (target is 12/06/2023). 

## Authors

Contributors names and contact info

* Joaquín de Navascués [@jdenavascues](https://twitter.com/jdenavascues)/[ORCID](https://orcid.org/0000-0002-5414-4056)
* Vinicius Dias Nirello [Google Scholar](https://scholar.google.com/citations?user=uMXPCs4AAAAJ))
* Aleix Puig-Barbé [@AleixPuig7](https://twitter.com/AleixPuig7)/[ORCID](https://orcid.org/0000-0001-6677-8489)

## Acknowledgments

Code snippets taken from Stack Overflow and other places are linked where they are used.
