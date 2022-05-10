# RMiXCR
The RMiXCR package is user friendly R interface to diminish the required bioinformatic skills for using the [MiXCR](https://mixcr.readthedocs.io/en/master/) software in the TCR recognition.
It is mainly focused in their use in the analysis of RNAseq bulk samples, facilitating its installation and usability by minimizing the inputs parameters, thus providing transparency and reproducibility of the results.

# Package Installation
```
install.packages("devtools")
devtools::install_github("elmerfer/RMiXCR")
```

# MiXCR installation
You can install RMiXCR in your favorite directory "the/path". After installation, you will see in your file structire the following file path directory "the/path/RMiXCR""
```
library(RMiXCR)
InstallMiXCR(where = "the/path")
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



