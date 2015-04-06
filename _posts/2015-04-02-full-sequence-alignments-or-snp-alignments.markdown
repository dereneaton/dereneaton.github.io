---
layout: post
title: "Should we discard invariable sites in phylogenetics?"
modified:
categories: phylogenetics
excerpt: "Computational speed is not a good excuse"
tags: [phylogeny, base-frequencies, invariable, raxml]
image:
  feature: header.png
  credit: Deren Eaton
date: 2015-04-02T20:12:50-04:00
comments: true
---

I've come across a number of studies in recent years in which a research group with perfectly good genomic sequence data interested in inferring phylogenetic relationships has elected to discard all invariable sites in the alignment and to instead use the raxml ascertainment bias setting which allows the analysis to be run with only SNP data. I don't recall ever reading a justification for this approach, but I always assumed that they had made the assumption that this method would be _faster_, since the data matrix being input to the program is much smaller. 

"In general I wouldn't use asc. bias. corr. if I had the invariable data, using the correction only makes sense if the invariable data has not been sampled/sequenced" -- Alexis Stamatakis (raxml google-group). 

And there is a memory consumption calculator on the Exelixis Lab web site as well that can be used to calculate how memory consuming your alignment will be. As Alexei says in the forum, the invariable sites are essentially a single column with weighted by their occurrence, so their memory consumption is miniscule. 

It seems that if memory is the limit on your analysis then you should be using examl instead of raxml, not using raxml on a subset of your data set. 

;)

:)

😙
👿
