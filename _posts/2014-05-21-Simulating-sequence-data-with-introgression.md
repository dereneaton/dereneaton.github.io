---
layout: post
title: "Simulating sequence data with introgression"
modified: 2014-05-21 16:58:29 -0400
category: phylogenetics
tags: [popgen, phylogenetics, simulation, introgression]
image:
  feature: header.jpg
  credit: Deren Eaton
  creditlink: 
comments: disqus
share: 
---


### Background

This script was developed to investigate how different tree shapes affect our ability to detect introgression between lineages. A tree topology with branch lengths is input as a newick string, and introgression events are designated with a string of values. It can simulate any number of loci of any given length, with no recombination within loci, but free recombination between them. The coalescent population parameter Ne is taken as an argument, and branch lengths are taken as Ne generations. The per-site mutation rate is fixed at 1e-9. The results are written to `stdout` in ".loci" format, which is similar to a multi-fasta format but has a line below each locus indicating the location of SNPs. It can be easily converted to alternative formats. The script has two dependencies, the Python packages [Egglib](http://egglib.sourceforge.net/) and [Numpy](http://www.numpy.org/). 

### Input arguments
The argument order is as follows:  

1. Number of loci
2. Locus length in bp
3. Effective population size (Ne)
4. Tree topology with branch lengths as a newick string 
5. Introgression events (see below)

The fifth argument, which is optional, lists introgression events as a string of five values, comma-separated, and inside of brackets. Multiple introgression events are separated by a "/" symbol. The five values in each string indicate:  

1. the donor lineage
2. the recipient lineage
3. time at which introgression starts
4. time at which introgression ends
5. probability (proportion) of migration (migrants)

Time is in coalescent units starting at 0.0 and increasing back in time. As an example, the string [C,D,0.20,0.25,1e-6] would indicate introgression from lineage D into C for 0.05 units of time. One coalescent time unit is equal to Ne generations, and the number of migrants per generation is 4 * Ne * m, where m is the last value passed in the introgression event argument. As a full example:  

### Calling the script
{% highlight python linenos %}
tree1="((((A:0.5,B:0.5):0.25,C:0.75):0.25,((D:0.5,E:0.5):0.25,F:0.75):0.25):2.0,((G:0.5,H:0.5):0.25,I:0.75):2.25);"
python simLoci.py 10000 100 1e6 '$tree1' [C,D,0.20,0.25,2.5e-7]
{% endhighlight %}

Or to enter multiple introgression scenarios:
{% highlight python linenos %}
python simLoci.py 10000 100 1e6 '$tree1' [C,D,0.20,0.25,2.5e-7]/[C,F,0.20,0.25,2.5e-7]
{% endhighlight %}

This would simulate 10K loci each 100bp in length on the input tree topology. It will allow 1 migrant per generation (4∗ 1e6 ∗ 2.5e−7=1) for 50K generations (0.05 ∗ 1e6=50K), for a total of 50K migrants into a population size of 1M.

{% highlight bash %}
>A         GAAGACCTATACGAATCTCCAAATATGGTGGATCACCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCGTG
>B         GAAGACCTATACGAATCTCCAATTATGGTGGATCACCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>C         GAAGACCTATACGAATCTCCAAATATGGTGGATCACCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>D         GAAGACCTATACGAATCTCCAAATATGGTGGATCACCACCGCCCCCCCAAAAGCGAAATGTCATCCACTACCCGACAGATGTAGGCTAAAGCTACTCCTG
>E         GAAGACCTATACGAATCTCCAAATATGGTGGATCACCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>F         GCAGACCTATACGAATCTCCAAATATGGTGGATCACCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>G         GAAGACCTATACGAATCTCCAAATATGGTGGATCAGCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>H         GAAGACCTATACGAATCTCCAAATATGGTGGATCAGCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
>I         GAAGACCTATACGAATCTCCAAATATGGTGGATCAGCACCGCCCCCCCAAAAGGGAAATGTCATCCACTACCCGACAGATGTAGGATAAAGCTACTCCTG
//          -                    -            *                 -                               -           -  |
>A         GTGCACCACATGCTATACTMATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGGGGTGTACCTCCATTTCGCC
>B         GTGCACCACATGCTATACTCATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGAMAGTCATCACGTGGGGTGTACCTCCATTTCGCC
>C         GTGCACCACATGCTATACTCATATTASAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACASTCATCACGTGGGGTRTACCTCCATTTCGCC
>D         GTGCACCACATGCTATACTCATATTASAACGGKCCGGGTCTGTTGCAATATATAGATATCCCACGGACASTCATCACGTGGGGTGTACCTCCATTTCGCC
>E         GTGCACCACATGCTATACTCATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGGGGTGTACCTCCATTTCGCC
>F         GTGCACCACATGCTATACTCATATTAGAACGGTCMGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGGGGTGTACCTCCATTTCGCC
>G         GTGCACCACATACTATACTCATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGTTGTGTACCTCCATTTCGCC
>H         GTGCACCACATGCTATACTCATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGTTGTGTACCTCCATTTCGCC
>I         GTGCACCACATGCTATACTCATATTAGAACGGTCCGGGTCTGTTGCAATATATAGATATCCCACGGACAGTCATCACGTGTTGTGTACCTCCATTTCGCC
//                    -       -      *     - -                                - *          **  -               |
>A         TGCCCGTAAACGCACATTTCACGACGTGTCGTATGGAAAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTC
>B         TGCCCGTAAACGCACATTTCACGACGTGTCGTATGGAAAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTC
>C         TGCSCGTAAACGCACATTTCACGACGTGTCGTATGGAAAAKGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTA
>D         TGCCCGTAAACGCACATTTCACGACGTGTCGTATGGAAAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAGCTCTGCAGCAAGCGTTCCTC
>E         TGCCCGTAAACGCACATTTCACGACGTGTCGTATGGAAAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTC
>F         TGCCCGTAAACGCACATTTCACGRCGTGTCGTATGGAAAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTC
>G         TGCCAGTAAACGCACATTTCACGACGTGGCGTATGGAGAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCWAACTCTCCAGCAAGCGTTCCTC
>H         TGCCAGTAAACGCACATTTCACGACGTGGCGTATGGAGAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTCCAGCAAGCGTTCCTC
>I         TGCCAGTAAACGCACATTTCACGACGTGTCGTATGCAGAAGGAACGGTGATAGCTGAAGCAACTGGTTTCGTATACCAAACTCTGCAGCAAGCGTTCCTC
//            -*                  -    *      - *  -                                    - -    *              -|
{% endhighlight %}