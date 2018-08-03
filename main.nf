import static groovy.json.JsonOutput.*
/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

def helpMessage() {
    log.info"""
    =============================================================
    csiro-crop-informatics/biokanga_align_paper  ~  version ${params.version}
    =============================================================
    Usage:

    nextflow run csiro-crop-informatics/biokanga_align_paper

    Default params:
    outdir      : ${params.outdir}
    publishmode : ${params.publishmode} [use 'copy' or 'move' if working across filesystems]

    Input params:
    Default settings and input specification suitable for test runs are read from conf/input.config
    Full pipeline run is triggered when executed with -params-file conf/input.json




    Trials/debugging params:
    trialLines  : ${params.trialLines}  [specify an int to subset input for trials, debugging]

    """.stripIndent()
}
// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//PARAMS
trialLines = params.trialLines

//ARRANGE INPUTS FOR PROCESSES
referencesLocal = Channel.create()
referencesRemote = Channel.create()
params.references.each {
  //Abbreviate Genus_species name to G_species
  it.species = (it.species =~ /^./)[0]+(it.species =~ /_.*$/)[0]
  //EXPECT TO HAVE SOME DATASETS WITH fasta
  if(it.containsKey("fasta")) {
    if((it.fasta).matches("^(https?|ftp)://.*\$")) {
      referencesRemote << it
    } else {
      referencesLocal << [it,file(it.fasta)]
    }
  }
}
referencesRemote.close()
referencesLocal.close()


process fetchRemoteReference {
  tag{meta.subMap(['species','version'])}
  label 'download'

  input:
    val(meta) from referencesRemote

  output:
    set val(meta), file("${basename}.fasta") into referencesRemoteFasta

  script:
    basename=getTagFromMeta(meta)
    //DECOMPRESS?
    cmd = (meta.fasta).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
    //TRIAL RUN? ONLY TAKE FIRST n LINES
    cmd += trialLines != null ? "| head -n ${trialLines}" : ""
    """
    curl ${meta.fasta} ${cmd} > ${basename}.fasta
    """
}

//Mix local and remote references then connect o multiple channels
referencesRemoteFasta.mix(referencesLocal).into{ references4rnfSimReads; references4kangaIndex; references4bwaIndex; references4bowtie2Index }


process indexReferences4rnfSimReads {
  tag{meta}
  label 'samtools'

  input:
    set val(meta), file(ref) from references4rnfSimReads

  output:
    set val(meta), file(ref), file('*.fai') into referencesWithIndex4rnfSimReads

  script:
  """
  samtools faidx ${ref}
  """
}

process rnfSimReads {
  tag{simmeta}
  label 'rnftools'

  input:
    set val(meta), file(ref), file(fai) from referencesWithIndex4rnfSimReads
    each nsimreads from params.simreads.nreads.toString().tokenize(",")*.toInteger()
    each length from params.simreads.length.toString().tokenize(",")*.toInteger()
    each simulator from params.simreads.simulator
    each mode from params.simreads.mode //PE, SE
    each distance from params.simreads.distance //PE only
    each distanceDev from params.simreads.distanceDev //PE only

  output:
    set val(simmeta), file("*.fq.gz") into reads4bwaAlign, reads4bowtie2align, reads4kangaAlign


  when:
    !(mode == "PE" && simulator == "CuReSim")

  script:
    tag=meta.species+"_"+meta.version+"_"+simulator
    simmeta = meta.subMap(['species','version'])+["simulator": simulator, "nreads":nsimreads, "mode": mode, "length": length ]
    len1 = length
    if(mode == "PE") {
      //FOR rnftools
      len2 = length
      tuple = 2
      dist="distance="+distance+","
      distDev= "distance_deviation="+distanceDev+","
      //FOR meta
      simmeta.dist = distance
      simmeta.distanceDev = distanceDev
    } else {
      len2 = 0
      tuple = 1
      dist=""
      distDev=""
    }
    """
    echo "import rnftools
    rnftools.mishmash.sample(\\"${tag}_reads\\",reads_in_tuple=${tuple})
    rnftools.mishmash.${simulator}(
            fasta=\\"${ref}\\",
            number_of_read_tuples=${nsimreads},
            ${dist}
            ${distDev}
            read_length_1=${len1},
            read_length_2=${len2}
    )
    include: rnftools.include()
    rule: input: rnftools.input()
    " > Snakefile
    snakemake \
    && for f in *.fq; do \
      paste - - - - < \${f} \
      | awk 'BEGIN{FS=OFS="\\t"};{gsub("[^ACGTUacgtu]","N",\$2); print}' \
      | tr '\\t' '\\n' \
      | gzip --stdout  --fast \
      > \${f}.gz \
      && rm \${f};
    done \
    && find . -type d -mindepth 2 | xargs rm -r
    """
}

process kangaIndex {
  label 'biokanga'
  label 'index'
  tag{dbmeta}

  input:
    set val(meta), file(ref) from references4kangaIndex

  output:
    set val(dbmeta), file(kangadb) into kangadbs

  script:
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    biokanga index \
    --threads ${task.cpus} \
    -i ${ref} \
    -o kangadb \
    --ref ${ref}
    """
}

process bwaIndex {
  label 'bwa'
  label 'index'
  tag{dbmeta}

  input:
    set val(meta), file(ref) from references4bwaIndex

  output:
    set val(dbmeta), file("${basename}*") into bwadbs

  script:
    dbmeta = ["species": meta.species, "version": meta.version]
    basename = getTagFromMeta(meta)
    """
    bwa index -a bwtsw -b 1000000000 -p ${basename}  ${ref}
    """
}

process bowtie2Index {
  label 'bowtie2'
  label 'index'
  tag{dbmeta}

  input:
    set val(meta), file(ref) from references4bowtie2Index

  output:
    set val(dbmeta), file("${basename}*") into bowtie2dbs

  script:
    basename=getTagFromMeta(meta)
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    bowtie2-build --threads ${task.cpus} ${ref} ${basename}
    """
}

process bwaAlign {
  label 'bwa'
  label 'samtools'
  label 'align'
  tag {alignmeta}

// no switch for soft clipping, but perhaps could set -L very high to disable:
// -L INT	Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query.
// If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied.

  input:
    set val(simmeta), file("?.fq.gz"), val(dbmeta), file('*') from reads4bwaAlign.combine(bwadbs) //cartesian product i.e. all input sets of reads vs all dbs

  output:
    set val(alignmeta), file('out.bam') into bwaBAMs

  when: //only align reads to the corresponding genome
    simmeta.species == dbmeta.species && simmeta.version == dbmeta.version

  script:
    dbBasename=getTagFromMeta(dbmeta)
    alignmeta = dbmeta + simmeta
    alignmeta.aligner = "bwa"
    if(simmeta.mode == 'SE') {
      """
      bwa mem -t ${task.cpus} ${dbBasename} 1.fq.gz | samtools view -bS > out.bam
      """
    } else {
      """
      bwa mem -t ${task.cpus} ${dbBasename} 1.fq.gz 2.fq.gz | samtools view -bS > out.bam
      """
    }

}

process kangaAlign {
  label 'biokanga'
  label 'align'
  tag {alignmeta}

 // params to explore: 2-3MM, indels, chimeric trimming


  input:
    set val(simmeta), file("?.fq.gz"), val(dbmeta), file(kangadb) from reads4kangaAlign.combine(kangadbs) //cartesian product i.e. all input sets of reads vs all dbs

  output:
    set val(alignmeta), file('out.bam') into kangaBAMs

  when: //only align reads to the corresponding genome!
    simmeta.species == dbmeta.species && simmeta.version == dbmeta.version

    //todo:EXPLORE PARAM SPACE
      // --microindellen 9 \
      // --minchimeric 50 \

  script:
    alignmeta = dbmeta + simmeta
    alignmeta.aligner = "biokanga"
    if(simmeta.mode == "SE") {
      """
      biokanga align \
      -i 1.fq.gz \
      --sfx ${kangadb} \
      --threads ${task.cpus} \
      -o out.bam
      """
    } else if(simmeta.mode == "PE"){
      """
      biokanga align \
      -i 1.fq.gz \
      -u 2.fq.gz \
      --sfx ${kangadb} \
      --threads ${task.cpus} \
      -o out.bam \
      --pemode 2
      """
    }
}

process bowtie2align {
  label 'bowtie2'
  label 'samtools'
  label 'align'
  tag {alignmeta}

  //explore parameter space: --local vs --end-to-end (default)

  input:
    set val(simmeta), file("?.fq.gz"), val(dbmeta), file('*') from reads4bowtie2align.combine(bowtie2dbs) //cartesian product i.e. all input sets of reads vs all dbs

  output:
    set val(alignmeta), file('out.bam') into bowtie2BAMs

  when:  //only align reads to the corresponding genome!
    simmeta.species == dbmeta.species && simmeta.version == dbmeta.version

  script:
    dbBasename=getTagFromMeta(dbmeta)
    alignmeta = dbmeta + simmeta
    alignmeta.aligner = "bowtie2"
    if(simmeta.mode == 'SE') {
      """
      bowtie2 -p ${task.cpus}  -x ${dbBasename} -U 1.fq.gz -p ${task.cpus} | samtools view -bS > out.bam
      """
    } else {
      """
      bowtie2 -p ${task.cpus} -x ${dbBasename} -1 1.fq.gz -2 2.fq.gz -p ${task.cpus} | samtools view -bS > out.bam
      """
    }
}

process rnfEvaluateBAM {
  label 'rnftools'
  tag{alignmeta}


  input:
    set val(alignmeta), file('out.bam') from bowtie2BAMs.mix(bwaBAMs).mix(kangaBAMs)

  output:
     set val(alignmeta), file(summary) into summaries

  script:
  // println prettyPrint(toJson(alignmeta))
  """
  rnftools sam2es -i out.bam -o - | awk -vOFS="\t" '\$1 !~ /^#/ {category[\$7]++};END{for(k in category) {print k,category[k]}}' > summary
  """

// rnftools sam2es OUTPUT header
// # RN:   read name
// # Q:    is mapped with quality
// # Chr:  chr id
// # D:    direction
// # L:    leftmost nucleotide
// # R:    rightmost nucleotide
// # Cat:  category of alignment assigned by LAVEnder
// #         M_i    i-th segment is correctly mapped
// #         m      segment should be unmapped but it is mapped
// #         w      segment is mapped to a wrong location
// #         U      segment is unmapped and should be unmapped
// #         u      segment is unmapped and should be mapped
// # Segs: number of segments
// #
// # RN    Q       Chr     D       L       R       Cat     Segs
}

process collateResults {
  label 'stats'
  executor 'local' //explicit to avoid a warning being prined. Either way must be local exec as no script block for this process just nextflow/groovy exec

  input:
    val collected from summaries.collect()

  output:
    file '*' into collatedResults

  exec:
  def outfileJSON = task.workDir.resolve('results.json')
  def outfileTSV = task.workDir.resolve('results.tsv')
  categories = ["M_1":"First segment is correctly mapped", "M_2":"Second segment is correctly mapped",
  "m":"segment should be unmapped but it is mapped", "w":"segment is mapped to a wrong location",
  "U":"segment is unmapped and should be unmapped", "u":"segment is unmapped and should be mapped"]
  entry = null
  entries = []
  i=0;
  TreeSet headersMeta = []
  TreeSet headersResults = []
  collected.each {
    if(i++ %2 == 0) {
      if(entry != null) {
        entries << entry
        entry.meta.each {k,v ->
          headersMeta << k
        }
      }
      entry = [:]
      entry.meta = it.clone()
    } else {
      entry.results = [:]
      it.eachLine { line ->
        (k, v) = line.split()
        //entry.results << [(k) : v ]
        entry.results << [(categories[(k)]) : v ]
        // headersResults << (k)
        headersResults << (categories[(k)])
      }
    }
  }
  entries << entry

  outfileJSON << prettyPrint(toJson(entries))

  //GENERATE TSV OUTPUT
  SEP="\t"
  outfileTSV << headersMeta.join(SEP)+SEP+headersResults.join(SEP)+"\n"
  entries.each { entry ->
    line = ""
    headersMeta.each { k ->
      line += line == "" ? (entry.meta[k]) : (SEP+entry.meta[k])
    }
    headersResults.each { k ->
      value = entry.results[k]
      line += value == null ? SEP+0 : SEP+value //NOT QUITE RIGHT, ok for 'w' not for 'u'
    }
    outfileTSV << line+"\n"
  }

}


process generatePlots {
  label 'rscript'
  label 'figures'

  input:
    file '*' from collatedResults

  output:
    file '*'

  script:
  '''
  #!/usr/bin/env Rscript

  #args <- commandArgs(TRUE)
  #location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
  if(!require(reshape2)){
    install.packages("reshape2")
    library(reshape2)
  }
  if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)
  }
  res<-read.table("results.tsv", header=TRUE, sep="\t");
  res2 <- melt(res, id.vars = c("aligner", "dist", "distanceDev", "mode", "nreads", "simulator", "species", "version","length"))
  pdf(file="results.pdf", width=16, height=9);
   ggplot(res2, aes(x=aligner, y=value,fill=variable)) +
   geom_bar(stat="identity",position = position_stack(reverse = TRUE)) +
   coord_flip() +
   theme(legend.position = "top") +
   facet_grid(simulator~mode~species);
  dev.off();
  '''
}


// /*
//   Generic method for merging read meta with db meta and aligner info
//  */
//  def getAlignMeta(meta, dbmeta) {
//   return (meta + dbmeta) //.aligner =
// }