---
layout: post
title: "Simulating raw RADseq data on a fixed tree"
modified: 2015-04-01
categories: radseq
excerpt: "fastq data for testing assembly"
tags: [radseq, phylogeny]
image:
  feature: header.png
  credit: Deren Eaton
date: 2015-03-31T18:24:18-04:00
comments: true
---

---------------------  

#### github link: ([_simrlls.py_](https://github.com/dereneaton/simrrls))
_This post has been modified for the updated v.1.04 script_.  

--------------------  

This script simulates raw (fastq format) sequence data
on a fixed species tree under the coalescent in a manner 
that emulates restriction-site associated DNA, with slight 
variations for different data types (e.g., RADseq, ddRAD, GBS, 
paired-end data). 

I originally developed this script to generate data for bug 
testing _pyRAD_ and to create example data sets for use in tutorials. 
The main purpose of this script is to test assembly methods
on raw data, if you are interested in simulating data for 
downstream analyses I would recommend instead using my 
([simLoci.py](/software/)) script, which simulates
_assembled_ RADseq-like data in a number
of usable formats (.phy, .nex, .geno, .migrate), 
and includes options to allow introgression between lineages.

Still, I figured a script to simulate raw data might be useful to 
others so I am making it available with some instruction.
To use the script you will need to first install the excellent [_Egglib_ 
Python package](http://egglib.sourceforge.net/), 
used here to perform the coalescent simulations, and it also requires
the common Python package _Numpy_.

__Calling the script__ -- _simrlls.py_ takes eight sequential arguments, 
each described in more detail below:

1.  Allow indels (decimal). The probability a mutation is a deletion.
2.  Allow mutations to restriction sites (int). 0=no locus dropout.
3.  Number of loci to sample (int; multiple of 100 or 1000)
4.  Number of individuals per tip taxon (int; minimum 1).
5.  Depth of sequencing (int, or (int,int)=(mean,std))
6.  Insert size (size selection window) (int,int)=(min,max)
7.  Data type (string; options = rad, ddrad, gbs, pairddrad, pairgbs)
8.  Prefix name for output files (string)

__Fixed arguments__--
For simplicity, 100 bp sequences are evolved on a fixed 12 taxon species 
tree. Similarly hard-coded is a sequencing error rate of 0.0005 mutations 
per site, and a mutation-scaled effective population 
size (theta) per tip taxon of 0.0014. The first restriction site
overhang is from PstI (TGCAG) and the second (used in ddRAD data) 
is from EcoRI (GAATT). These fixed parameters can be easily 
changed in the first few lines of the script.
The outgroup taxon in the tree "X" is not included in the output 
data, but is used in the simulation to polarize mutations relative 
to an outgroup for creating indels and missing data. 

![simtreeimage](/images/setupsims.png)

__arg 1__ -- Point mutations relative to the outgroup taxon "X"
will be converted to a deletion with the probability set here (e.g., 0.01). 
If 0, no indels are present in the data set. 

__arg 2__ -- This option allows locus dropout to occur with a rate that is 
scaled relative to (theta * max size selection window * length of restriction
recognition sites). If the the maximum length of size selected fragments
is 500, then the locus is dropped if a mutation occurs relative to the 
outgroup such that either of the restriction recognition sites occurs 
within the length of this fragment. The locus is also dropped if 
a mutation occurs within the restriction site(s). Locus dropout can 
be turned off (set to 0). 

__arg 3__ -- Nloci is the number of loci that will be generated for the
first sampled individual in the outgroup taxon after excluding loci that
would have experienced locus dropout in this sample. If locus dropout is 
turned off then all other samples will have the same number of loci. 
If locus dropout is turned on, all locus dropout
occurs by mutations relative to this sample. 

__arg 4__ -- Ninds is the number of sampled individuals from each tip 
taxon with divergence among individuals scaled by theta, which
 is fixed at 0.0014. 

__arg 5__ -- Sequencing depth for diploid individuals. If an int is 
entered then it is rounded up to the nearest even number and the 
two alleles if present are sampled with equal frequency. E.g., depth=20
would mean each allele is sampled 10 times. Alternatively two 
comma-separated int values can be entered, upon which each _allele_ will 
be sampled from a normal distribution with values _mean,stdev_. An 
example with values 10,2 could be that at one locus allele\_1 is sampled
11 times and allele\_2 is sampled 9 times. 

__arg 6__ -- Insert size. 
This option will affect the rate at which locus dropout occurs
and can also be used to introduce merging of paired reads, or overlap 
of short GBS reads if the minimum insert size is <0. 
As long as the lower range of the insert size is >0 no overlap 
of paired end reads or GBS reads occurs. 
Insert sizes are uniform randomly selected from within the 
user supplied size window. The total size selected fragments are 
200bp + insert size. 

__arg 7__ -- This arguments determines the data type. 
Depending on the choice the script will put cut sites and barcodes in 
the proper ends of sequences so that they emulate RAD, GBS, or ddRAD 
data, single or paired-end. The choice also effects how locus dropout occurs.

#### Calling the script
{% highlight bash %}
python simrrls.py 0.01 1 1000 10 20 400,800 rad simrads
{% endhighlight %}

_The following will be printed to the screen_  
{% highlight bash %}  
simulating rad data
1000 loci at 20X coverage
10 samples per taxon across 12 tip taxa
indels arise at frequency of 0.0100 per mutation
locus dropout = True
size selection window = 400,800
sequencing error rate = 0.0005
theta=4Nu= 0.0014
{% endhighlight %}

__Output__ -- Two data files are created, for our example these
are `simrads_R1_.fastq.gz` and `simrads.barcodes`. The first is a
compressed fastq file with sequence data and the 
latter is a text file mapping 
barcodes to sample names. If you select a paired-end data 
type it will create two sequence files, one with "\_R1\_" 
in the name and the other with "\_R2\_". 
With the data and a barcode map the data can then be analyzed
in either _pyRAD_ or a similar program such as _stacks_. 

_simrads.barcodes (only the first few lines)_  
{% highlight bash %}
1A0     CATCAT
1A1     TTTTGG
1B0     GGAGTA
1B1     TAAGTT
1C0     TGAGGT
1C1     TTAAGT
1D0     AATAGG
1D1     TGGGGG
2E0     TAAAAG
2E1     TGATTA
{% endhighlight %}

_simrads\_R1\_.fastq.gz (only the first few lines)_  
{% highlight bash %}
@lane1_locus1_0_R1_0 1:N:0:
ATAGGGTGCAGTTTGGGAAAATTGAGTGAACCCCGGCGTGATCTAGCCCGCGCAGAAGAGGATGACCGCTGCGCCTCTCACTCACTTATGCTACTAATAC
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@lane1_locus1_0_R1_1 1:N:0:
ATAGGGTGCAGTTTGGGAAAATTGAGTGAACCCCGGCGTGATCTAGCCCGCGCAGAAGAGGATGACCGCTGCGCCTCTCACTCACTTATGCTACTAATAC
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@lane1_locus1_0_R1_2 1:N:0:
ATAGGGTGCAGTTTGGGAAAATTGAGTGAACCCCGGCGTGATCTAGCCCGCGCAGAAGAGGATGACCGCTGCGCCTCTCACTCACTTATGCTACTAATAC
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@lane1_locus1_0_R1_3 1:N:0:
ATAGGGTGCAGTTTGGGAAAATTGAGTGAACCCCGGCGTGATCTAGCCCGCGCAGAAGAGGATGACCGCTGCGCCTCTCACTCACTTATGCTACTAATAC
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
{% endhighlight %}


#### Some more examples: 
For an ideal data set used in testing _pyRAD_ I would simulate data with equal coverage across loci, without indels, and with no locus dropout:
{% highlight bash %}
python simrlls.py 0.00 0 1000 1 20 100,400 rad simrads_testing
{% endhighlight %}

If I wanted to emulate real data I would input variable coverage across loci, include indel variation and locus dropout:
{% highlight bash %}
python simrlls.py 0.01 1 1000 5 10,5 100,400 rad simrads_realistic
{% endhighlight %}

And if I wanted to test the ability of _pyRAD_ (or pyRAD in combination with PEAR read merging) to filter out overlapping paired-end reads I would simulate paired-end data with a size selection window smaller than the paired read lengths:
{% highlight bash %}
python simrlls.py 0.00 0 1000 1 20 -50,200 pairddrad pairddrad_w_merging
{% endhighlight %}
