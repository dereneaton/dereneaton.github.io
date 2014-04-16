---
layout: post
title: "New static Jekyll website"
description: "What is Jekyll and why use it?"
category: web
tags: [jekyll,web]
image:
  feature: header.jpg
  credit: 
  creditlink: 
comments: true
share: 
---

### What is jekyll and why use it?
Welcome to www.dereneaton.com, my new static Jekyll website. [Jekyll](http://www.jekyllrb.com) is a Ruby parsing engine for converting plain text into static webpages, meaning you can you can create web pages or blog posts using markdown, a simple language for formatting text, code, figures and links. Here I am using the [_Minimal Mistakes_](http://mademistakes.com/articles/minimal-mistakes-jekyll-theme/) Jekyll theme for the overall page design. Over time I will tweak and customize it to better fit my needs, with all changes saved to my github repository. Anyone can fork the site at any time to begin creating their own Jekyll site based on my build. 

And this is perhaps the best feature of a Jekyll site, GitHub offers a free web hosting service powered by Jekyll, meaning you can host your entire site on Github and make changes to it using version control. I considered alternative static page generators, particularly those written in Python like Pelican and Nikola, but even a cursory glance at the downloads or example pages makes clear that, at least in terms of theme development and refinement, Jekyll is miles ahead. 

### Creating this website
I followed the instructions here ([_MinimalMistakes_](http://mademistakes.com/articles/minimal-mistakes-jekyll-theme/)) to create my site, which was as easy as cloning a git repository, editing a few markdown files in a text editor, and pushing the changes to a git repository. I was then able to transfer the host name of my old WordPress site to this one following the instructions here ([link](https://help.github.com/articles/my-custom-domain-isn-t-working)).

### Transferring a domain name
The only tricky part of this was finding the IP address to enter into the domain hosting site (I use namecheap.com), this is found with the 'dig' command. Your site on github will be created with your github username and ending with github.io, so mine was dereneaton.github.io. To transfer this to dereneaton.com I simply had to create a file in my repository call CNAME and entere "dereneaton.com" into it, then point the hosting service to the github repository. This is described in the link above. The IP address you will enter is the one after github.map.fastly.net:

{% highlight bash %} 
dig dereneaton.github.io +nostats +nocomments +nocmd

;www.dereneaton.github.io.	  IN	      A
www.dereneaton.github.io. 1549	  IN	      CNAME	github.map.fastly.net.
github.map.fastly.net.	  30	  IN	      A		199.27.76.133

{% endhighlight %}

In terms of creating download links, posting blog entries, and generally designing the site in a way I find pleasing, I find Jekyll to be extraordinarily more intuitive than WordPress. For the non-GUI users out there, comfortable with a terminal and plaintext, trust me, creating static webpages is the way to go. There are many other advantages too, as you can find with a simple google search, including increased security and speed.  Hopefully this post is useful to others interested in setting up a simple static blog/website. 

