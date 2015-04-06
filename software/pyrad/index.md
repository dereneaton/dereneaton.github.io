---
layout: page
title: "pyrad"
date: 2015-3-31
modified:
excerpt: 
tags: [software, pyrad]
share: false
image: 
  feature: header.png
  credit: Deren Eaton
---

### Description

The benefit of _pyRAD_ over most alternative methods for analyzing RADseq-like data comes in its use of an alignment-clustering method (_vsearch_) that allows for the inclusion of indel variation which improves identification of homology across highly divergent samples. For this reason _pyRAD_ is commonly employed for RADseq studies at deeper phylogenetic scales, however, it works equally well at shallow scales.  

### Flexibility

_pyRAD_ is intended for use with any type of restriction-site associated DNA. It currently supports RAD, ddRAD, PE-ddRAD, GBS, PE-GBS, EzRAD, PE-EzRAD, 2B-RAD, nextRAD, and can be extended to other types. Below is the download link as well as a number of tutorials. There is a “general use” tutorial which explains how to install and setup input files, and also tutorials with example data sets for different data types and analyses. The software was initially described in the following publication ([preprint](http://biorxiv.org/content/early/2013/12/03/001081), [journal](http://bioinformatics.oxfordjournals.org/content/early/2014/03/20/bioinformatics.btu121 ))  

### Downloads (now on github)  
You can download a stable release version, or if you are comfortable with _git_ you may clone the repository. Whichever you choose, please check for updates frequently, as I generally make bug fixes or updates weekly (see [changelog](https://github.com/dereneaton/pyrad/commits/master)).  

+  [Current stable release](https://github.com/dereneaton/pyrad/releases)
+  [git Development version](https://github.com/dereneaton/pyrad/)

#### General Use Tutorial
+  [Full tutorial v.3.0](http://nbviewer.ipython.org/gist/dereneaton/af9548ea0e94bff99aa0/pyRAD_v.3.0.ipynb)  

#### Recommended example Tutorial  
+  [SE RAD Tutorial v.3.0](http://nbviewer.ipython.org/gist/dereneaton/1f661bfb205b644086cc/tutorial_RAD_3.0.ipynb)  

#### Advanced tutorials  
+  [SE ddRAD Tutorial v.3.0.4](http://nbviewer.ipython.org/dc6241083c912519064e/tutorial_ddRAD_3.0.4.ipynb)
+  [PE ddRAD Tutorial v.3.0.4](http://nbviewer.ipython.org/dc6241083c912519064e/tutorial_pairddRAD_3.0.4.ipynb)
+  [PE ddRAD w/ merged reads Tutorial v.3.0.4](http://nbviewer.ipython.org/gist/dc6241083c912519064e/tutorial_pairddRAD_3.0.4-merged.ipynb)
+  [SE GBS (or ezRAD) Tutorial v.3.0.4](http://nbviewer.ipython.org/gist/dereneaton/9d12ff5ab6584c5ceafa/tutorial_GBS_3.0.ipynb)  
+  [PE GBS (or ezRAD) w/ merged reads Tutorial v.3.0.4](http://nbviewer.ipython.org/gist/dereneaton/1f661bfb205b644086cc/PE-GBS_empirical.ipynb)  
+  D-statistics tutorial()  

(For now, contact the google group forum for additional supported datatypes: 2B-RAD, nextRAD, etc.)  

--------------------------  

#### Google group questions forum ([link](https://groups.google.com/forum/#!forum/pyrad-users))

---------------------------  

#### [Previous version tutorials link](http://nbviewer.ipython.org/gist/dereneaton/af9548ea0e94bff99aa0)  


