---
layout: post
title: "Simulating RADseq data on a fixed tree"
description: "a script using Egglib"
modified: 2014-04-15 17:41:56 -0400
category: phylogenetics
tags: [radseq,phylogeny,popgen]
image:
  feature: header.jpg
  credit: Deren Eaton
  creditlink: 
comments: true
share:
---


---------------------  

### Download the script
+ [simRRLs.py](/downloads/simRRLs.py)

--------------------  

### Background
I originally developed this script to simulate data for bug testing _pyRAD_,
but more recently I've used several variants of it to simulate data 
with different amounts of divergence between samples, different proportions of 
missing data, and with or without gene flow between lineages. 
I figured it may be useful to others so I am making 
it available here, with some instruction.

The script sets up a species tree and simulates sequence data on this tree under a 
coalescent model, then dresses up the data to act as if it had been generated as a
reduced representation genomic library (RRL), such as RADseq, ddRAD or GBS. 
To use the script you will need to first install the excellent [_Egglib_ 
Python package](http://egglib.sourceforge.net/), 
used here to perform the coalescent simulations, and it also requires
the common Python package _Numpy_.

### Calling the script
_simRLLs.py_ takes seven arguments, some of which are further described below:

1.  Allow indels (decimal). Probability a mutation is a deletion.
2.  Allow mutations to restriction sites (1/0).
3.  Number of loci to sample (int; multiple of 100 or 1000)
4.  Number of individuals per tip taxon (int; minimum 1).
5.  Size selection of fragments (int,int)
6.  Data type (string; options: rad, ddrad, gbs, pairddrad, pairgbs)
7.  Prefix name for output files (string)

#### Sequences are evolved on this 12 taxon species tree
Ideally one would be able to pass as an argument any topology on which to evolve the sequences,
but for now there is simply one topology hard-coded into the script. 
This can easily be changed by editing the script. The outgroup taxon "X" 
is not included in the output data, but is used in the simulation 
to polarize mutations relative to an outgroup. 

![simtreeimage](/images/setupsims.png)

##### arg 1  
Any point mutations relative to the outgroup taxon "X" will be converted to a deletion
with the probability set here. If 0, no indels are present in the data set. 

##### arg 2  
This option allows mutations to arise in a region (such as the restriction site)
that will cause locus dropout. The size of this region, and thus the relative rate at
which locus dropout occurs, is set as an integer. Any mutations arising in this region
relative to the outgroup will cause reads to be excluded from the library. 
This type of "locus dropout" is commonly considered a problem with empirical 
RADseq data. By simulating sequences this way they 
will have missing data that is phylogenetically correlated similar to empirical data sets. 

##### arg 5 
This is mainly used for developing filtering/trimming options in _pyRAD_. If the fragment
length is very small relative to the sequence length then sequences will include 
Illumina adapters, and some reads (GBS or paired-end) can include the barcode sequences
where two reads overshoot each other in the reverse direction. As long as the lower range
is higher than 200 this has no effect on the results. Fragment lengths are randomly selected from within the size window. 

##### arg 6
Locus dropout will differ between RAD and the other methods since RAD data have 
a restriction cut site on only one side of sequences, whereas the other library
preparation methods have restriction sites on both ends of sequences, and thus greater
probability (assuming cut sites occur with equal frequency) of locus dropout. This script does not accurately reflect the true vagaries of how these methods will yield missing data. This option will format the data with restriction sites in the proper places on the sequences, and with single or paired end data. The proportion of missing data is controlled by arg 2 alone. 

#### Calling the script
{% highlight python linenos %}
python simRRLs.py 0 0.1 1000 10 400,800 rad simrads
{% endhighlight %}

#### Output
Two data files are created, for our example these are `simrads_R1.fastq.gz` and `simrads.barcodes`. The first is a compressed fastq file with sequence data and the latter is a text file mapping 
barcodes to sample names. If you select a paired-end data type it will create two sequence 
files, one with "R1" in the name and the other with "R2". With the data and a barcode map
the data can then be analyzed in either _pyRAD_ or a similar program such as _stacks_. 


#### Adding introgression between lineages
Adding a line of code in the script like below can introduce gene flow
between different lineages of the tree. In this case we will
set that migration occurs between time=0.5 and time=1.0
from taxon number 3 (named "1C0") into taxon number 2 (named "1B0"), 
with a frequency of 4\*N\*0.0001 migrants per generation. 
Remember time is counted backwards from the present (coalescent),
so migration is zero at time=0, then is starts at time=0.5,
and we set gene flow back to zero at time=1.0. 

{% highlight python linenos %}
paramSet.changePairwiseMigrationRate(0.5, 2, 3, 4*N*0.0001)
paramSet.changePairwiseMigrationRate(1.0, 2, 3, 0.0)
{% endhighlight %}

