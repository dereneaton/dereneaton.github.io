---
layout: page
title: "Software"
date: 2015-3-31
modified:
excerpt:
tags: [software]
share: false
image: 
  feature: header.png
  credit: Deren Eaton
---

### [_pyRAD_](/software/pyrad)  
Software to assemble _de novo_ RADseq loci from restriction-site associated sequence data (RAD,ddRAD,GBS,PEddRAD,PEGBS). Alignment clustering allows for indel variation to better align highly divergent samples, merge and trim methods canbe employed to rescue overlapping paired end data, and reverse complement clustering is available to improve GBS assemblies of short overlapping contigs. 

### [_simRRLs_](/software/simRRLs/)  
Python script to simulate fastQ formatted RAD, ddRAD or GBS data on a fixed species tree under a coalescent model. Generating this type of data is useful for exploring the capabilities of _pyRAD_, _STACKS_, or related programs for assembling RADseq data sets with varying amounts of divergence between lineages, and with indel variation. 

### [_simLoci_](/software/simLoci/)  
Python script to simulate sequence data on an input topology under a coalescent model with arguments to allow migration between lineages. In contrast to _simRRLs.py_, above, this script outputs aligned loci that do not require assembling. This is useful for investigating the affects of introgression on phylogenetic inference, or for testing metrics of detecting introgression between lineages. 

