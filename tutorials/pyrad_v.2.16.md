---
layout: page
permalink: /tutorial/pyrad_v.2.16/
title: "Tutorial pyRAD_v.2.16"
tags:
    - python
    - notebook
    - pyrad
---


### Deren A. R. Eaton  

#### Contact: deren.eaton@yale.edu  

#### Help forum: [google group link](https://groups.google.com/forum/#!forum/pyrad-users)

---------------

## Contents:

#### [1. Why use _pyRAD_?](-1.-Why-use-_pyRAD?-)  ####

#### [2. Installation](#2.-Installation)  

#### [3. Input files](#3.-Input-files)  

#### [4. Running _pyRAD_](#4.-Running-pyRAD)  

#### [5. RADseq/RRL datatypes](#5.-RADseq/RRL-datatypes)  

#### [6. Example Analyses](#6.-Example-Analyses)  

#### [7. Output data formats](#7.-Output-data-formats)  

#### [8. Hierarchical clustering](#8.-Hierarchical-clustering)  

#### 9. Parallel execution, and how to improve speed


----------------


# 1. Why use _pyRAD_? #

Reduced-representation genomic sequence data (e.g., RADseq, GBS, ddRAD) are
commonly used to study population-level research questions and consequently most
software packages for analyzing RADseq data are designed for data with little
variation across samples. Phylogenetic analyses typically include species with
deeper divergence times (more variable loci across samples) and thus a different
approach to clustering and identifying orthologs will perform better.

_pyRAD_, the software pipeline described here, is intended to maximize
phylogenetic information across disparate samples from RAD, ddRAD, or GBS data,
hereafter referred to generally as RADseq. Written in Python, the program is
flexible and easy to use. The code is human-readable, open-source and may be
modified to meet users specific needs. With respect to constructing phylogenetic
data sets, the largest difference between _pyRAD_ and alternative programs such
as _Stacks_, is in the way that loci are clustered. _Stacks_ clusters reads by
grouping those with less than N differences between them, but does not allow for
indels, such that the presence of even a single indel will cause a frameshift
resulting in many base differences between otherwise highly similar sequences.
_pyRAD_ instead uses a measure of sequence similarity based on a global
alignment of sequences through the program USEARCH. Alignment clustering allows
for indels, and should perform better in data sets including more distantly
related samples.

_pyRAD_ was written with the aim that it be easy to use, fast, and flexible. It
is easy in the sense that installation does not require any difficult libraries
or extensions, and that it requires few arguments to the command line. It is
fast by offerring a number of heuristics to speed clustering of very large data
sets, or long (paired-end) reads. And it is flexible in that it has built-in
methods to account for different library types (RAD, ddRAD, GBS, paired-ddRAD,
paired-GBS) as well as quality filtering, including overlap among paired or
reverse-complemented reads, and the presence of adapter sequences. Thus, data do
not need to be pre-screened or filtered before analysis in _pyRAD_.

Two fast clustering methods: hierarchical-clustering for large data sets, and
split-clustering for paired-end ddRAD data are described further below. Using
these methods population or phylogenetic scale data sets with even many hundreds
of individuals can be assembled in relatively short time by distributing across
multiple processing cores. Finally, _pyRAD_ can be used to perform D-statistic
tests on RADseq loci to test for historical introgression.

An open access pre-print describing _pyRAD_ and comparing it with _Stacks_ is
available here: [link](http://biorxiv.org/content/early/2013/12/03/001081)

## 2. Installation

_pyRAD_ is available for Linux and Mac (but not Windows). It can be downloaded
at the following [link](http://dereneaton.com/software). A few Python packages 
and external software are required as dependencies. I list
how you can find and install each below:

+ python 2.5+
+ python-numpy
+ python-scipy
+ muscle
+ usearch v.7.0.1+ (32 or 64 bit versions)

--------------------


`Numpy` and `Scipy` are common Python libraries required for numerical and
scientific calculations. On __Ubuntu Linux__ these can be installed using 
apt-get commands:

{% highlight bash %}
  sudo apt-get install python-numpy
  sudo apt-get install python-scipy
{% endhighlight %}

On a __Mac__ there are a few ways to install these, depending on the version of
your OS. I suggest you search 'Mac install scipy', or Scipysuperpack, or Conda,
for installing these two packages.

To check that they installed properly open Python in a command
terminal by typing "python" and enter the import arguments below. 
If everything is installed properly you will get no error messages.

{% highlight python %}
  ## in Python
  import numpy
  import scipy
  import scipy.stats
{% endhighlight %}

---------------

The other two __required__ pieces of software are `MUSCLE` and `USEARCH`, both
available at [http://www.drive5.com](www.drive5.com). Once downloaded you can move 
the executables into your `$PATH` so they can be called from the command line, or, 
if you don't know what that means, simply take note of the location where you
store the executable files, as you will need to pass this information to _pyRAD_.

Sometimes the downloaded executable of `USEARCH` does not have the correct
permissions set for execution (I don't know why). If you get an error when you
execute the program saying you do not have permission to use it, then run the
following unix command to change permissions of the file:

{% highlight bash %}
  chmod 555 usearch_v.7.0.1090_i86linux32
{% endhighlight %}
----------

__Important note__: All versions of _pyRAD_ v.1.70 and above are compatible with
`USEARCH v.7.0.1090` and above, not with older versions.

-----------


## 3. Input files

### The barcodes file

The barcodes file is a simple table linking barcodes to samples. 
Each line should have one name and one barcode, separated by a
tab (important: not spaces).

{% highlight bash %}
 sample1     ACAGG
 sample2     ATTCA
 sample3     CGGCATA
 sample4     AAGAACA
{% endhighlight %}

Barcodes can be of varying lengths. If your data are already de-multiplexed
(sorted) by samples into separate files then a barcode file is not needed.

### The params file 

A template parameter file, typically named __params.txt__, should be created using the
(__-n__) option in _pyRAD_, like below. This file lists all of the options or
paramater settings necessary for an analysis. Each line contains a parameter value 
(a number or character string) followed by
__any number of spaces or tabs and then two hash marks (“##”)__, after which the
parameter is described, or comments can be added. In parentheses is listed the step
of a _pyRAD_ analysis affected by the parameter. If a line is left blank the
default value is used. I highly recommend beginning all analyses by creating an
empty template params.txt file with the following command:

{% highlight python %}
  ~/pyrad_v.2.1/pyRAD -n
{% endhighlight %}

-------------  

This creates a file with the default settings which you will want to edit
with a text editor.

--------------  

The following section describes the params file line by line. The first 14 lines
are required, all lines following this are optional and used for more advanced
analyses or data filtering. For each input line I list several possible
formats:

### required lines 1-14 ###

+ __Line 1__: The working directory for your analysis. Three different examples
are shown below. This is where all of the output files from the analyses will be
placed.

{% highlight bash %}
                           ## 1. uses current as working directory  
        ./                 ## 1. uses current as working directory  
        /home/user/RAD/    ## 1. uses set location as working directory  

{% endhighlight %}

--------------------

+ __Line 2__: The location of the raw fastq formatted Illumina sequence data.
Use the wildcard character \* to select multiple files in the same directory.
Paired-end data files should be in pairs that differ only by the presence of a
"_R1_" and "_R2_" in their names. Data files can be gzipped.

{% highlight bash %}
        ./raw/*.fastq       ## 2. path to raw data files
        raw/*.fastq.gz      ## 2. path to raw data files (gzipped)
{% endhighlight %}

---------------------------------

+ __Line 3__: The location of the barcodes file (described above). If your data
are already de-multiplexed then lines 2 and 3 can be left blank, but line 18
(see below) must be entered. Example:

{% highlight bash %}
        ./barcodes.txt        ## 3. path to barcodes file
        mydata.barcodes       ## 3. path to barcodes file
{% endhighlight %}  

-----------------------

+ __Line 4__: The command to call USEARCH. If USEARCH is in your $PATH, then
simply list the executable file name here. Otherwise, enter the full path to
USEARCH that is required for you to call it from the command line. Examples:

{% highlight bash %}
         /user/bin/usearch7.0.1090_i86linux32   ## 4. full path
         usearch7.0.1090_i86linux32             ## 4. in system PATH
{% endhighlight %}

----------------------------

+ __Line 5__: The command to call _muscle_. Similar to above.

{% highlight bash %}
        muscle                       ## 5. path to call muscle
        muscle3.8.31_i86linux32      ## 5. path to muscle long name
{% endhighlight %}


-------------------------------

+ __Line 6__: The sequence of the restriction recognition site overhang. If your
data are ddRAD list the common cutter (attached to first reads) followed by the
rare cutter (attached to second reads), separated by a comma. If you're not sure
of the sequence of your cutters look at your raw data files for the common
sequence found at the left side of your reads.  Cutters can contain ambiguous
bases (Example: CWGC). Examples:

{% highlight bash %}
         TGCAG                  ## 6. PstI cutsite overhang C|TGCAG
         CWGC,AATT              ## 6. ddRAD example for apeKI, EcoRI
{% endhighlight %}

-------------------------------

+ __Line 7__: The number of processors (cores) to use. Entering a value that is
more than the number of available processors will yield slower results. Parallel
processing in _pyRAD_ works by executing individual samples separately on
different processors, so you cannot gain any speed by setting a number greater
than the number of samples in your data set. Example:

{% highlight bash %}
         8                      ## 7. N processors
{% endhighlight %}

-------------------

+ __Line 8__: The minimum depth of coverage to make a statistical base call at
each site in a cluster. Depending on the inferred error rate, five is typically
the minimum depth at which a base can be called. You may want to use a higher
minimum depth, like 10. If your data are very low coverage you can instead call
majority consensus base calls on reads with coverage as low as 2 by setting
option 31 (see below).

{% highlight bash %}
         10                     ## 8. mindepth for base calls
{% endhighlight %}

----------------------

+ __Line 9__: The maximum number of low quality, undetermined (“N”) sites in
step 2 filtered sequences. This number should be selected with respect to the
clustering threshold. A low clustering threshold (.85) can allow more
undetermined sites; whereas a high clustering threshold (.95) should exclude
reads with too many undetermined sites or they will affect clustering.

{% highlight bash %}
         4                      ## 9. maxN: max number of Ns in reads
{% endhighlight %}

--------------------------------

+ __Line 10__: The similarity threshold to use for alignment clustering in
USEARCH, entered as a decimal. This value is used for both within-sample
clustering and across-sample clustering. I (more or less) require you to use the
same threshold for both, and recommend doing so. Given the way the clustering
algorithm works it only makes sense (to me) to do so, otherwise reads that do
not cluster together within a sample will later cluster together when you do
across-sample clustering, which is undesirable. During clustering the cutsite
bases remain attached to the reads, thus all loci will have an invariable set of
4-6 bases on one side. For 100 bp reads containing a 5-bp cutter and with have a
5 bp barcode trimmed off this results in 95 bp of data, 5 of which are
invariable. Thus, if you want to allow 5 base differences between reads (90/95 =
.947) you should use a setting of .94. Example:

{% highlight bash %}
        .90             ## 10. clustering threshold
{% endhighlight %}

------------------------------

+ __Line 11__: The type of Reduced representation library. There are five
options, currently. Entered as all lowercase lettering: rad, ddrad, gbs,
pairddrad, pairgbs. The different treatments of these data types are described
below in section RADseq/RRL datatypes.

{% highlight bash %}
        rad               ## 11. Datatype
        pairddrad         ## 11. Datatype
{% endhighlight %}

---------------------------------

+ __Line 12__: Minimum taxon coverage. The minimum number of (ingroup) samples
with data for a given locus to be retained in the final data set. If you enter a
number equal to the full number of samples in your data set then it will return
only loci that have data across all samples. If you enter a lower value, like 4,
it will return a more sparse matrix, including any loci for which at least four
samples contain data. Example:

{% highlight bash %}
        4               ## 12. minCov
{% endhighlight %}

----------------------------

+ __Line 13__: Maximum number (or proportion) of shared polymorphic sites in a
locus. Enter a number, or a decimal with the prefix “p” (e.g., p.10 for 10%).
This option is used to detect potential paralogs, as a shared heterozygous site
across many samples likely represents clustering of paralogs with a fixed
difference rather than a true heterozygous site. Example:

{% highlight bash %}
        4               ## 13. maxSharedH.
{% endhighlight %}

------------------------------

+ __Line 14__: Prefix name for final output files. I usually enter a name
indicative of the other parameters from the lines above. Example:

{% highlight bash %}
        c90d10m4p4      ## 14. output name prefix
{% endhighlight %}

-------------------------------------

### optional lines 15-35 ###

These options are not necessary for all analyses and can be left empty in the
params file, in which case _pyRAD_ uses their default values. They mostly
pertain to options for additional quality filtering, for paired-end analyses, or
for increasing the speed of analyses.

--------------------------------


+ __Line 15__: Subset selector for steps 2-7. If, for example, you have data for
a number of individuals from two different clades, and for one clade all of the
sample names begin with the letter “A”, then you can select only these samples
to analyze by entering the shared prefix "A” on this line.

{% highlight bash %}
        A          ## 15. selects prefix subset of files
{% endhighlight %}

-----------------------------------

+ __Line 16__: Add-on (outgroup) taxa selector. The MinCov parameter on line 12
tells _pyRAD_ to keep only loci with data from at least n ingroup samples. If
you wanted to maximize coverage within a certain clade you could set the other
taxa as "outgroups", and only loci with at least MinCov ingroup samples will be
kept, and any additional matching of ’outgroup’ taxa to those will be included
in the data, but not count toward the minimum number of samples. List ’outgroup’
taxa names here separated by commas.

{% highlight bash %}
        outg1,outg2            ## 16. outgroup selector for step 7
{% endhighlight %}

-----------------------------------

+ __Line 17__: Exclude taxa. When you are constructing data sets at step 7 you
can exclude certain taxa in order to maximize coverage among other included
samples. For example you may want to exclude taxa that are found to have very
few data. List sample names comma separated here. Example:

{% highlight bash %}
        sample3,sample4        ## 17. Exclude from step 7
{% endhighlight %}

------------------------------------

+ __Line 18__: If your data are already de-multiplexed into separate files for
each barcode then you can skip step 1 in the analysis and go straight to step 2.
Now, you must enter the location of your sorted fastq files on this line. The
reads should have the barcode trimmed off. If they also have the cutsite trimmed
off then enter the ’@’ symbol at the beginning of the line to indicate this.
Files can be compressed (.gz) or not. Examples:

{% highlight bash %}
        ~/RAD/fastq/*.fastq    ## 18. Sorted files
        @~/RAD/fastq/*.fastq   ## 18. sorted and stripped files
{% endhighlight %}

------------------------------------

+ __Line 19__: The maximum number of mismatches allowed in a barcode during 
de-multiplexing. Default is 1, I don’t generally recommend going above this.

{% highlight bash %}
        1                    ## 19. max barcode mismatches
{% endhighlight %}

-------------------------------------

+ __Line 20__: Phred Quality score offset (usually 33, but sometimes 64). Base
calls with a phred score below 20 are converted to Ns. If you want this value to
be more strict you can change the offset to be higher. For example, to raise the
minimum score to 30 set the offset to 43. If left blank the default is 33.

{% highlight bash %}
        33                   ## 20. min Quality score step 2
{% endhighlight %}

-------------------------------------

+ __Line 21__: Strictness of filtering in step 2. (0) means no filtering for
barcodes, adapters and cut sites, only correct for quality scores of base calls.
(1) looks for cutsite+adapters and trims them out. (2) tries to detect these
while allowing some errors in them, thus enforcing the most strict filtering.
For most data sets that do not include many overlapping and short fragments
options 0 or 1 is recommended.

{% highlight bash %}
        0                    ## 21. strictness of Q filter
{% endhighlight %}

-------------------------------------

+ __Line 22__: _a priori_ error rate and heterozygosity. Step 4 of the analysis
jointly estimates the error rate and heterozygosity, and step 5 uses these
values to make base calls. If for some reason you wish to skip the error rate
estimate (e.g., your samples are polyploid), then you can enter an _a priori_
estimate of the error rate here. Caution: you should generally let the program
estimate it by using step 4, as the error rate here is being measured on data
that have already passed through some filters, so error rates estimated
elsewhere will not be correct. Across several data sets I typically find error
rates between 0.0001 - 0.002, depending on length and sequence chemistry.
Example:

{% highlight bash %}
        0.0005,0.001        ## 22. E and H, respectively
{% endhighlight %}

--------------------------------------

+ __Line 23__: maximum number of “N”s in a consensus sequence. Clustering across
samples can be affected by the number of “N”s, and so this number should be
chosen with respect to line 6 and the lengths of reads. Default is 4.

{% highlight bash %}
        5                  ## 23. MaxN in consensus seqs
{% endhighlight %}

---------------------------------------

+ __Line 24__: maximum number Hs (heterozygous sites) in a consensus sequence.
Among several methods to exclude potential paralogs or highly repetitive genomic
regions, you can set a maximum number of heterozygous sites in a consensus
sequence. Default is 4.

{% highlight bash %}
        5                  ## 24. MaxH in a consensus seq
{% endhighlight %}

-----------------------------------------

+ __Line 25__: Allow only 2 haplotypes in a consensus sequence. In diploid
organisms, after correcting for sequencing errors, there should be at most two
haplotypes making up any consensus genotype call. by default this option will
exclude loci with more than 2 alleles. If set to 0 it will not apply a filter.
More advanced polyploid options will be implemented in the future. Example:

{% highlight bash %}
        1                  ## 25. diploid filter max haplos=2
{% endhighlight %}

-----------------------------------------

+ __Line 26__: Maximum number of SNPs in a final locus. This can remove
potential effects of poor alignments in repetitive regions in a final data set
by excluding loci with more than N snps in them. For paired data you can set max
SNPs in the first read and second read separately, or use one value for the
pair, Examples:

{% highlight bash %}
        10               ## 26. max SNPs in a final locus
        8,12             ## 26. max SNPs in 1st and 2nd pair read
{% endhighlight %}

-----------------------------------------

+ __Line 27__: Maximum number of insertions/deletions. If only one value,
applies to within-sample clusters. If two values, first is within-sample
clusters, second is across-sample clusters. For paired-data, enter four values:
withinclusterread1, withinclusterread2, acrossclustersread1,
acrossclustersread2. Defaults are 3,6,99,99 respectively. Examples:

{% highlight bash %}
        3                ## 27. max indels in a within-sample cluster
        3,10             ## 27. max within-sample, max for across-sample cluster
        3,6,10,10        ## 27. max
{% endhighlight %}

within_read1,within_read2,across_read1,across_read2

-------------------------------------------

+ __Line 28__: random number seed. Randomization is used to sort the input order
of sequences before clustering in steps 3 and 6. This can have an effect on
results. If you set a random number seed then analyses should be repeatable.

{% highlight bash %}
        112233           ## 28. random number seed
{% endhighlight %}

-------------------------------------------


+ __Line 29__: Allow overhanging ends of reads in final data set. If reads are
different lengths or overlap to different degrees. Options, (x,x) effect the left and
right ends of a single-end sequence or of a paired contig. Alternatively, if you 
enter four values (x,x,x,x) this effects the left and right of the first read, followed
by the left and right of the second read, respectively. 0 means no trimming, 1 means 
trim to the site where at least four samples have data (not "-"), and 2 means trim to the
site where all samples at that locus have data. 

{% highlight bash %}
        1,1              ## 29. trim to informative sites
        0,0              ## 29. no trim 
        0,2,2,0          ## 29. trim far ends of paired reads
        2,2,2,2          ## 29. trim both ends of paired reads
{% endhighlight %}

--------------------------------------------

+ __Line 30__: Output formats (u,s,a,n). u = unlinked snps, s = snps, a =
alleles, n = nexus. See output/formatting section.

{% highlight bash %}
        a,n              ## 30. outputs formats
{% endhighlight %}

-------------------------------------------

+ __Line 31__: Call simple consensus (majority base) on clusters with depth less
than mindepth, (allows you to set option 8 as low as 2).

{% highlight bash %}
        1                ## 31. call majority base on low depth sites
{% endhighlight %}

--------------------------------------------------

+ __Line 32__: Check for overlap of paired reads. Enter minimum overlap of
paired reads. Overlapping reads are merged and written to a separate file, and a
file is created with only non-merged reads.

{% highlight bash %}
        30               ## 32. minimum overlap (nbases) of pairs
{% endhighlight %}

------------------------------------------------------

+ __Line 33__: Keep trimmed sequences. This option is only for very messy single
end GBS or ddRAD data, where a size selection method did not work properly such
that many fragments were shorter than the read length and thus the adapters were
frequently sequenced. Instead of discarding sequences with adapters that are
trimmed this option includes those that are still longer than X.

{% highlight bash %}
         0               ## 33. don't keep trimmed sequences
         50              ## 33. keep trimmed longer than 50
{% endhighlight %}

--------------------------------------------------------

+ __Line 34__: maximum depth of a within-sample cluster. The default (empty) is
max(meandepth+2*SD, 500), meaning the higher value of either the mean plus two
times the standard deviation of cluster depth, or 500. If instead you want to
set an absolute value enter it here. Or enter a ridiculously high value for no
maxdepth filtering.

{% highlight bash %}
                           ## 34. default max depth option
         50000             ## 34. a set max depth
         999999999999999   ## 34. no max depth
{% endhighlight %}

---------------------------------------------------------

+ __line 35__: minimum number of copies of de-replicated reads used for
clustering. Reads are dereplicated, or collapsed, into a single read that
records its number of occurrences before clustering. The default is to not
exclude any data. But for exploratory analyses this can be useful because the
speed of within-sample clustering can be greatly improved by excluding singleton
reads, or reads that occurred less than N times. In general I do not recommend
discarding data for your final analysis, but for data sets with deep coverage it
has little effect on results and greatly improves clustering times.

{% highlight bash %}
         1                 ## 35. exclude singletons
         5                 ## 35. exclude reads occurring less than 5 times
{% endhighlight %}

----------------------------------------------------------------

+ __Remaining lines__: Hierarchical clustering. If performing hierarchical (two-
step) clustering the assignment of individuals to groups and the minimum cluster
size within groups are listed here. Each group is listed on a separate line, it
is OK to leave a line blank. Each line will contain three elements, each
separated by a space: (1) a group name (e.g., 'A' or '1'), these names are
arbitrary, just choose something different for each group. (2) the minimum size
of a cluster within this group, any cluster smaller than this will not be
clustered across-groups. (3) the names of samples/individuals in this group.
Sample names can be individually listed comma separated, or selected using a
wildcard selector (e.g., \*). Example:

{% highlight bash %}
        A 4 1A0,1B0,1C0,1D0
        B 4 2*
        C 4 3*
{% endhighlight %}


__________________________

## 4. Running _pyRAD_

_pyRAD_ should be called from the directory in which it comes (i.e., with its
other dependency .py files). For example, to call _pyRAD_ if it had been
downloaded to a user’s home/ directory, type the command below, which will print
the help menu.

{% highlight python %}
  ~/pyrad_v.2.1/pyRAD -h
{% endhighlight %}

-------  

There are five main options to the command line interface of _pyRAD_:

+ -h (help screen)
+ -n (generate new params file)
+ -p (designate params file)
+ -s (choose steps)
+ -d (D-test input file)


------------------


The parameter input file (params.txt) should be created using the -n option and
then filled in based on the type of analysis you plan to do. If the params file
is properly filled out then the entire analysis – converting raw RADseq data
into a assembled de novo loci – can be done by simply entering:

{% highlight python %}
   ~/pyrad_v.2.1/pyRAD -p params.txt
{% endhighlight %}

This will perform all steps (1-7) of the analysis, as described below. If,
however, you wish to perform only one, or a few step at a time, then you can
select these using the -s (steps) option. To select only the first step:

{% highlight python %}
  ~/pyrad_v.2.1/pyRAD -p params.txt -s 1
{% endhighlight %}

And to run steps 2-7 you could enter:

{% highlight python %}
  ~/pyrad_v.2.1/pyRAD -p params.txt -s 234567 
{% endhighlight %}

After finishing one assembly, you could rerun step 7 with the same or a
different params.txt file indicating different step 7 filtering options.

{% highlight python %}
  ~/pyrad_v.2.1/pyRAD -p params.txt -s 7
{% endhighlight %}

When first learning new software it often takes a while to learn how to properly
fill in the correct parameter settings and options. This can be espeicially
frustrating with computationally intensive analyses, where one might have to
wait many hours (days) before learning that they made a simple mistake. The -s
option is an attempt to reduce this annoyance by breaking up the analysis into
seven discrete jobs, each of which is performed sequentially. In this way, you
will know whether step 1 worked before moving on to step 2. And if step 2 fails,
you will retain the results of step 1, and can start again from where the error
arose. I also suggest running the example data sets to make yourself familiar
with the steps of an analysis and what the output files and statistics look like
before beginning your own analysis.

### The seven steps described

+ __Step 1__: de-multiplexing. This step uses information from the barcodes file
to separate sequences from your raw fastq files into a separate file for each
sample. These are placed in a new directory within your working directory called
“fastq/”. For paired-end data it is necessary that the raw data file names
follow a specific format: The first read files must contain "\_R1\_" in the
name, and the second read files must be identical to the first read files but
with "\_R2\_" in place of "\_R1\_". Here is an example pair of input files:

{% highlight bash %}
  yourfilename_R1_001.fastq
  yourfilename_R2_001.fastq
{% endhighlight %}

_A reminder_, if your data are already demultiplexed this step can be skipped. You
will not need to enter a location for your raw data in the params file, but
instead you must enter the location of your demultiplexed files into the proper
location on line 18 of the params file. 

+ __Step 2__: filtering. This step uses the Phred quality score recorded in the
FASTQ data files to filter low quality base calls. Sites with a score below a
set value are changed into “N”s, and loci with more than the number of allowed
“N”s are discarded. Files are written to the “edits/” directory with the suffix
“.edit”. It also implements a number of optional filters. For paired-end data
you can separate out overlapping (merged) sequences (line 32), and you can filter
for Illumina adapter sequences (line 21).

+ __Step 3__: within-sample clustering. This step first dereplicates the filtered
sequences from step 2, recording the number of times each unique read is
observed. These are then clustered using USEARCH, recording if they match within
a set sequence similarity. Sequences that clustered together are then aligned
and written to a new file in the “clust.xx/” directory with the ending
“.clustS.gz”. This is typically the most time consuming step of an analysis. 

+ __Step 4__: error-rate and heterozygosity estimates. This step uses ML
to jointly estimate error rate and heterozygosity from
the base counts in each site across all clusters. 
Results are written to the Pi_estimate.txt file in the stats/ directory.

+ __Step 5__: create consensus sequences. Using the mean error rate and
heterozygosity estimated in step 4, this step creates consensus sequences for
each cluster. Those which have less than the minimum coverage, more than the
maximum number of undetermined sites, or more than the maximum number of
heterozygous sites, or more than the allowed number of alleles, are discarded.
If two alleles are present the phase of heterozygous sites are retained for the
consensus sequences.

+ __Step 6__. Consensus sequences are clustered across samples using the same
settings as in step 3. If heterozygous, one allele is randomly sampled and used
in this step, although both alleles are retained in the final data set.

+ __Step 7__. Alignment, detection of paralogs, output of human readable fasta-
like file (.loci), and concatenated loci (.phy), other optional formats can be
specified (.nex, .snps, .alleles), and final statistics are written to .stats
file. This step is relatively fast, and can be repeated with different values
for options 12,13,16,17 to create different data sets that are optimized in
coverage for different samples.

## 5. RADseq/RRL datatypes

The following reduced representation library data types are supported by
_pyRAD_, and I list beside each my general experience with analyzing empirical
data sets. My quality judgements are based on the frequency with which users
have generated data sets that include many short fragment reads with adapter
sequences, or which require extra types of analyses (e.g., reverse complement
matching for overlapping or paired-end GBS.) This may change as library
preparations of GBS and ddRAD libraries are improving. I describe each data type
and how _pyRAD_ treats them differently below.

      Type                Quality                      Per-sample speed
      __________________________________________________________________________
      RAD                 (very good)                     Fast (1-3 days)
      nextRAD             (very good)                     Fast (1-3 days)
      ddRAD               (usually good)                  Fast (1-3 days)
      paired ddRAD        (sometimes difficult/messy)     Average (1-7 days)
      GBS                 (sometimes difficult/messy)     Slow (3-14 days)
      paired GBS          (often difficult/messy)         Very slow (5-21 days)

![svg]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_57_0.svg)


##### datatypes

In general single-end Illumina data that employ a good size selection window
will produce very good and usable data by any of the library prep methods above.
Problems tend to arise when the size selection is too small compared to the sequencing
length, and especially when this coincides with paired-end sequencing.

Poor size selection has little effect on traditional RAD or on single end ddRAD,
as you can see above, except to include adapter sequences in the reads which can
be filtered out. It has a more significant effect on paired ddRAD, and a very
significant effect on single or paired GBS data. This is because short fragments
will be sequenced from both ends of these types of data, leading to sequence
overlap.

In paired ddRAD the first and second reads come in pairs, and _pyRAD_ can test
for overlap of these reads using the merge overlap option. This will remove
merged reads and write them to a separate file. These can be analyzed separately
from the non-overlapping reads if they make up a large proportion of your data
set. Recommendations for this type of analysis are made in the paired ddRAD
tutorial.

For GBS data the problem is more difficult. In single end data the Illumina
adapter can ligate to either end of fragments, meaning short single-end
fragments can have overlap leading to duplication in the data set. _pyRAD_ can
correct for this by performing reverse-complement clustering. This is much
slower than normal clustering but can still yield high quality data sets. 
Recommendations for this type of analysis are made in the GBS tutorial.

For paired-end GBS data the merge overlap option can merge overlapping reads,
however the designation of a read as first versus second is completely 
interchangeable, so reverse complement clustering is also required in this case,
and because reads are longer it usually takes very long.

_________________________________  

## 6. Example Analyses 

The links below are to a number of analyses including example data sets that you
can download and try for yourself.

##### Basics:

1. [Simple RADseq analysis](http://goo.gl/WT5Fhc)

2. [Simple paired ddRAD analysis](http://goo.gl/gB8Zq9)

3. [Simple GBS analysis]

4. [Simple ddRAD analysis]


##### Advanced Examples:

5. [paired ddRAD with overlap and Hierarchical clustering of many
individuals](http://nbviewer.ipython.org/gist/dereneaton/c18bff4ba8bf532ec14b)

6. [GBS analysis with short fragment rescue (poor size selection)](http://nbviewer.ipython.org/gist/dereneaton/9d12ff5ab6584c5ceafa)  

7. [D-statistic tests for introgression]

______________________________________  


## 7. Output data formats

By default _pyRAD_ return the final data in two formats: individual loci
(__.loci__ file) and concatenated loci (__.phy__ file), in addition to outputting a
file with the loci that did not pass the final step of paralog filtering
(__.excluded_loci__), which you can use to fine tune your final parameter values
used in step 7. 

There are a number of additional formats possible as well. These include nexus
(__.nex__), haplotypes/alleles (__.alleles__), snps (__.snps__), unlinked
snps (__.unlinked\_snps__), and structure format (__.str__). These additional formats 
will be output if indicated on line 30 of the params file by a shorthand letter,
i.e, a for alleles, n for nexus, u for unlinked snps, s for snps, k for structure. 

### *.loci ###

This is a custom format that I use to output the individual loci and indicate
variable sites. For each locus it provides the sequence data for each sample and
indicates the location of variable sites or snps. Custom scripts can easily
parse this file for loci containing certain amounts of taxon coverage or
variable sites. Also it is the most easily readable file for assuring that your
analyses are working properly. Example below:

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_65_0.png)

For paired-end data the two linked loci are shown separated by a 'xxxx'
separator.

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_67_0.png)

### *.phy ###

This is a phylip formatted data file which contains all of the loci from the
.loci file concatenated into a supermatrix, with missing data for any sample
filled in with N's. This format is used in RAxML among other phylogenetic
programs.

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_69_0.png)

### *.nex ###

This is a nexus formatted data file which contains all of the loci from the
.loci file concatenated into a supermatrix, but printed in an interleaved
format, with missing data for any sample filled in with N's, and with data
information appended to the beginning. This format is used in BEAST among other
phylogenetic programs.

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_71_0.png)

### *.alleles ###

This is the same as the .loci format but instead of ambiguity codes for the
consensus sequences the two alleles are printed out for each sample at each
locus. Output of polyploid (>2) alleles are not yet supported.

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_73_0.png)

### *.snps ###

This format lists only the variable sites in each locus with a space between
SNPs from different loci. An underscore represents an invariable locus. Loci are
in the same order as in the .loci file. Paired loci are treated as a single
locus, meaning SNPs from the two loci are not separated in this file (linked).

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_75_0.png)

### *.unlinked_snps ###

Each column contains one SNP sampled from one locus. If multiple SNPs in a
locus, SNP sites that contain the least missing data across taxa are sampled, if
equal amounts of missing data, they are randomly sampled.

![png]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_77_0.png)

### *.excluded_loci ###

The excluded loci file is in the same format as the .loci file but shows the
loci that were filtered out in step 7. The filter is listed below each locus. %P
indicates the locus was filtered as a paralog because the samples shared a
heterozygous site over more than the allowed number of samples. %D indicates the
locus had more than one sequence from the same individual. %S indicates there
were more SNPs detected in the locus than the set maximum allowed number. %I
means there were more indels in the locus than the max allowed number.

{% highlight python %}
>A     CCRACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCACCGCCCCCGCGATCT
>B     CCRACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCACCGCCCCCGCGATCT
>C     CCRACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCACCGCCCCCGCGATCT
>D     CCRACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCACCGCMCCCGCGATCT
>E     CCRACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCAGCACCCCCGCGATCT
>F     CCAACACAGAGGGCCAATAGACACGTAGCTCGAGATCTCAGCAACCCCGCGATCT
//%P     *                                       *--          |
>A     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGACCGAGAACTAGGT
>B     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGACCGAGAACTAGGT
>C     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGACCGAGAACTAGGT
>D     CGTACAATAGGATGCACCGTTCGTCAAGTTGTTGACGGTGGACCGAGAACTAGGT
>D     CGTACAATAGGATGCACCGTTCGTCAAGTTGTTGACGGTGGACCGAGAACTAGGT
//%D                           -     -                        |
>A     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGACCGAGAACTAGGT
>B     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGGACGGACGGACGGA
>C     CGTACAATAGGATGCACCGTTCGTAAAGTTCTTGACGGTGGGACGGACGGACGGA
>D     CGTACAATAGGATGCACCGTTCGTCAAGTTGTTGACGGTGGACCGAGAACTAGGT
//%S                           -     -          **  *******  *|
...
{% endhighlight %}

______________________________________  

## 8. Hierarchical clustering  

Clustering across samples (step 6 of the analysis) can be quite fast for data
sets with few sampled individuals, but for studies with many (>100) samples this
step can take very very long (weeks), depending on the total number of unique
sequences. Hierarchical clustering splits this job up into multiple smaller jobs
that each run on a separate processing core, allowing clustering of even very 
large data sets with hundreds of individuals in as little as a few days. 

The greatest speed increase comes from excluding loci at the first step of
clustering that have significant missing data within a subgroup. By excluding
singletons, or loci with little shared data early in the clustering process time
is not spent comparing them to all of the other more distantly related samples.
This is particularly relevant if you plan to only use loci in your final data
set that recover data shared across many individuals.

![svg]({{ site.baseurl}}notebooks/pyrad_v.2.1_files/pyrad_v.2.1_82_0.svg)

The figure above illustrates the general concept. Individuals are grouped _a
priori_ into a small number of clades or populations, in this case __a__, __b__,
and __c__. The assignment of individuals does not have to be perfectly correct
in terms of knowing their true phylogeny, rather this is an approximation or
heuristic, where those grouped together are expected to share more data than
those which are not grouped together.

In this case each clade has four sampled taxa. To group our samples into these
three clades we assign them to groups in the bottom lines of the params file.
Each group is listed on a separate line, it is OK to leave blank lines. Each
line will contain three elements, each separated by a space: (1) a group name
(e.g., 'A' or '1'), these names are arbitrary, just choose something different
for each group. (2) the minimum size for a cluster within this group, any cluster
smaller than this will not be clustered across-groups. (3) the names of
samples/individuals in this group. Sample names can be individually listed comma
separated, or selected using a wildcard selector (e.g., \*).

The sample names being matched to are the prefix names in the
"__.consens.gz__" files located in the clust directory.  For example, imagine we
have the samples below.

{% highlight bash %}
sample1.consens.gz
sample2.consens.gz
sample3.consens.gz
sample4.consens.gz
sample5.consens.gz
sample6.consens.gz
sample7.consens.gz
sample8.consens.gz
sample9.consens.gz
{% endhighlight %}

Imagine also that we know samples 1-4 are individuals of species A, samples 5-8
are individuals of species B, and sample9 is from a more distant outgroup
species. We could cluster them using the following settings:

{% highlight bash %}
    A 4 sample1,sample2,sample3,sample4
    B 2 sample5,sample6,sample7,sample8
    C 1 sample9
{% endhighlight %}

The minimum cluster size option tells _pyRAD_ to only cluster loci across clades
that have at least N matches within a clade. In this way loci that are
singletons, or shared by very few taxa are not clustered across clades.  In the
case of the outgroup sample C, we entered 1 for the minimum size, and thus all
reads from this sample are clustered to the other two clades. In contrast, we
chose 4 for clade A, and 2 for clade B. This means only loci that have data for
all four samples in clade A will be clustered with clades B and C, and only loci
with data from at least 2 individuals in clade B will be clustered with clades A
and C.

Following this logic, you can assign individuals to clades in a way that speeds
clustering but also maximizes the amount of data retained in your data set. If a
clade has very few individuals in it, or those individuals have less data
overall than other samples then you may want to set a low minimum cluster size
to discard as little data as possible. In general, though, you only get dramatic
speed increases if at least some clades have minimum cluster size set >1.

There is an example tutorial implementing hierarchical clustering and which
provides more detail in the Examples section above.

___________________________________

## 9. Parallel execution, and how to improve speed 

There are a number of ways to make _pyRAD_ analyses run faster than the 
default settings, but there is usually a trade-off in the amount of data retained 
as these heuristics utilize ways of excluding data with very low coverage within or 
across samples. 

1. __Distribute jobs across many processors__ on multiple computers or computing nodes. 
_pyRAD_ executes steps 2-5 on each individual sample separately, and it is only in 
step 6 that they must be combined to cluster _across_ samples. Thus, the earlier steps, 
most importantly step 3 (within-sample clustering), can be performed 
simulataneously on as many different computing processors as you have available. This 
is easiest to do using the subset selector option in the params file. One thing to note, 
however, is that step 5 should not be executed until all samples have completed step 4, 
so that the mean estimate of error rate and heterozygosity across all samples can be 
used to make consensus base calls. Once step 4 is complete, step 5 can be run across
multiple processors. 
As an example, say you have a computing cluster with 4 nodes that each have 8 processors. 
You can submit a job to each node that selects 8 samples using the subset selector, 
thus running 32 jobs simultaneously. 

2. __Use the minderep option__. When this option is set above 1 (the default) 
it excludes dereplicated reads that occurred less than 'minderep' times. Because 
a single sequencing error, or uncalled base (N) can cause a sequence to uniquely
dereplicate, singletons usually represent a _very large_ proportion of the reads
that are compared to other reads during clustering. Clustering is ordered in _pyRAD_
with reads that dereplicate the highest being first, and thus more likely to be seeds, 
and singletons coming last. Excluding singletons will not effect how high copy sequences
cluster to one another, but will decrease the number of low copy clusters that reach 
sufficient depth to make statistical base calls and thus be retained as consensus 
sequences. Using this option excludes greater amounts of data with longer reads, 
because there are typically more unique reads. I recommend using at least minderep=2
for _exploratory_ analyses because it greatly improves the speed, but for full 
analyses you should always use the default minderep=1 if you have time to run 
it on your computing resources. 

3. __Use hierarchical clustering__. As described above, for any reasonably large data set
this can yield huge increases in speed during step 6 (across-sample clustering), 
with minimal effects on missing data. 

____________________________________  

#### Running _pyRAD_ on a computing cluster

+ remember to load the necessary modules (python, numpy & scipy).
+ _pyRAD_ does not use MPI but rather runs parallel jobs sequentially.
+ keep an eye on memory usage. Each process of _pyRAD_ does not use 
a lot of memory (at most 2GB, but usually much less). Only in 
step 4, the most memory intensive step, should you worry about using
fewer processors at once, since 50 processes simultaneously
may exceed your memory limit. 
+ distribute jobs across separate nodes by running them separately 
and using the subset selector. 

