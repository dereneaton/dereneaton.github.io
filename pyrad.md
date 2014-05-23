---
layout: page
permalink: /software/pyrad/
title: "pyRAD"
modified: 2014-03-22 23:11
tags: [pyrad,software]
image:
  feature: header.jpg
  credit: Deren Eaton
  creditlink: 
share: 
---

---------------------   

#### Description

The benefit of _pyRAD_ over some alternative methods for analyzing RADseq-like data comes in its use of an alignment-clustering method (_usearch_) that allows for the inclusion of indel variation which improves identification of homology across highly divergent samples. For this reason _pyRAD_ is commonly employed for RADseq studies at deeper phylogenetic scales, however, it works equally well at shallow scales.  

#### Flexibility

_pyRAD_ can analyze RAD, ddRAD, GBS, paired-end ddRAD and paired-end GBS data sets. It also has funtionality for trimming/removing adapter sequences and merging overlapping paired end reads. A two-step clustering method makes clustering even very large data sets with hundreds of samples feasible and fast, and a “split clustering” method for paired-end ddRAD data dramatically speeds assembly of paired loci.  

#### Continued development
Below you can find the current release, as well as older versions. There is a “general use” tutorial which explains how to install and setup input files, as well as specific tutorials with example data for different data types and analyses. The software is described in the following publication ([preprint](http://biorxiv.org/content/early/2013/12/03/001081), [journal](http://bioinformatics.oxfordjournals.org/content/early/2014/03/20/bioinformatics.btu121 ))  

## Current stable release download:  

+  [pyRAD_v.2.11.zip](/downloads/pyrad_v.2.11.zip)
+  [Tutorial\_pyRAD\_v.2.01+](/tutorial/pyrad_v.2.1/)

### Example data set notebooks (v.2.0+):  
+  [RADseq tutorial](http://nbviewer.ipython.org/gist/dereneaton/1f661bfb205b644086cc)  
+  [paired-end ddRAD tutorial](http://nbviewer.ipython.org/gist/dereneaton/c18bff4ba8bf532ec14b)  
+  [GBS tutorial](http://nbviewer.ipython.org/gist/dereneaton/9d12ff5ab6584c5ceafa)
+  ddRAD tutorial

+  D-statistic tests tutorial

--------------------------  

#### Google group questions forum ([link](https://groups.google.com/forum/#!forum/pyrad-users))

---------------------------  

#### git development ([version](https://code.google.com/p/pyrad/)) and change log ([link](https://code.google.com/p/pyrad/source/list))

__previous stable releases:__    
[pyrad_v.2.10.zip](/downloads/pyrad_v.2.10.zip)  
[pyrad_v.2.01 tutorial](http://nbviewer.ipython.org/gist/dereneaton/af9548ea0e94bff99aa0)  
[pyrad_v.2.01.zip](/downloads/pyrad_v.2.01.zip)  
[pyrad_v.2.00.zip](/downloads/pyrad_v.2.0.zip)  
[pyrad_v.1.64.zip](/downloads/pyrad_v.1.64.zip)  
[pyrad\_tutorial\_v.1.64.pdf](/downloads/pyrad_v.1.64.pdf)  
[pyrad_v.1.4.zip](/downloads/pyrad_v.1.4.zip)  
