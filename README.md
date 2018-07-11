# Introduction

* Expand on bread wheat genome challenges
    * polyploid
    * repeats

## Existing aligners 

how detailed do we need to be?
different alignment approaches. cover pseudo alignment? ie kallisto, salmon? 


# Methods

* Reference genome/s (probably one is enough?)
    * Zimin/Salzberg
    * Clavijo
    * one would be enough, but in a pipeline there should be none or very little overhead of having 2
    * can we quantify 'complexity'? ie repeats/duplications/syntenic regions with limited polymorphisms... ideally from existing work)
        * Worth doing even though this can be only indicative - really long (almost identical) repeats will be missing from the assemblies! 
        * would any of Stuart's hamming stats be of use? 
        * k-mer frequency analysis either on the assembly or some raw data (either CS or RAC8754 PE+MP) - should be quick, with KMC can go up to read length (max k=255). When a plot illustrating that includes various species, it can be very informative
        * [genome mappability score](https://academic.oup.com/bioinformatics/article/28/16/2097/323484) or [similar](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377)
* Align: (what parameter space? can't be too complete)
    * DNA reads (bioplatforms? - selected highest coverage cultivars?)
    * simulated reads
      * can we use biokanga simreads?
      * and/or some third-party?
      * Parameter space to be explored: error rates should be close to reported empirical 
    * alignment parameters e.g. max substitutions
      * could go slightly outside the range used when aligning reads from/to wheat 
      * soft clipping / chimeric trimming and their influence on speed, align rate, false positive rate...
      * PE or SE? If PE do we need to explore the various PE modes, how much can they contribute?
* Measure:
    * Standard parameters (time/cores blah blah)
    * %aligned in different classes (unaligned, unique etc)
    * Accuracy for simulated reads 
    * SNP calling? Typically this is done with a different tool to the aligner - are we introducing SNP calling as well? unsure. could work well with bioplatforms data.



## Simulated reads procedure

Much of this can be taken care of by using [https://github.com/karel-brinda/rnftools](https://github.com/karel-brinda/rnftools)

```
for each assembly #may be informative to show arabidopsis, rice, barley, wild emmer wheat, bread wheat 
  for each read_length #fix?
    for each err rate 
      simulate n-fold coverage reads (proportional to assembly/genome size) 
      for each param settings in {default, ..., ...} 
        align reads to ref using each tool
        quantify speed, accuracy
``` 



# Results

## Simulated reads

I can imagine a figure showing the fate of reads under different aligners. A sort of snakey diagram showing where they end up, and whether it's the best place or not. 

## bioplatforms reads

Summary of different SNPs reported (ie what is the overlap in the SNPs identified using the different aligners - would need to use the same SNPs caller on the BAM output I guess)? I don't actually know what the differences would look like. Would be nice to have an example or two of where a difference is. Maybe this is really down to parameters (ie how many reads you need to trust the call?)


# Discussion


