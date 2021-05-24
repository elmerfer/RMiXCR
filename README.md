# RMiXCR
The RMiXCR package is user friendly R interface to diminish the required bioinformaic skils for using the [MiXCR](https://mixcr.readthedocs.io/en/master/) software in the TCR recognition.
It is mainly focused in their use in the analysis of RNAseq bulk samples, facilitating its installation and usability by minimizing the imputs parameters, thus providing transparency and reproducibility of the results.

# Package Installation
```
install.packages("devtools")
devtools::install_github("elmerfer/RMiXCR")
```

# MiXCR installation
```
library(RMiXCR)
InstallMiXCR()
```

# Capabilities
It will provide
* easy installation of MiXCR software
* Transparent RNAseq analysis for TCR from V, D, J genes rearrangements.
* Couple with R data structures to facilitate its analysis
* Improved export capabilities to Excel(r) files
* Analitical capabilities in estimation diversity, gene usage an several new statistical capabilities

# Developers
* Elmer A. Fernàndez. (PhD), main author and idea. CIDIE-UCC-CONICET.
* Guadalupe Nibeyro. developer
* Veronica Baronetto developer, and analitical tools.



