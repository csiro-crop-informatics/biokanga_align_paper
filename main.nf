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
    Please modify/specify in conf/input.config

    nreads  : ${params.simreads.nreads}  

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

process rnfSimReads {
  tag{simmeta}
  label 'rnftools'

  input:
    set val(meta), file(ref) from references4rnfSimReads
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
  tag{dbmeta}

  input:
    set val(meta), file(ref) from references4bwaIndex

  output:
    set val(dbmeta), file("${ref}*") into bwadbs

  script:
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    bwa index -a bwtsw ${ref}
    """
}

process bowtie2Index {
  label 'bowtie2'
  tag{dbmeta}

  input:
    set val(meta), file(ref) from references4bowtie2Index

  output:
    set val(dbmeta), file("${basename}*") into bowtie2dbs

  script:
    basename=getTagFromMeta(meta)
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    bowtie2-build ${ref} ${basename}
    """
}

process bwaAlign {
  label 'bwa'
  label 'samtools'
  label 'align'
  tag {alignmeta}

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
      bwa mem -t ${task.cpus} ${dbBasename}.fasta 1.fq.gz | samtools view -bS > out.bam
      """
    } else {
      """
      bwa mem -t ${task.cpus} ${dbBasename}.fasta 1.fq.gz 2.fq.gz | samtools view -bS > out.bam
      """
    }

}

process kangaAlign {
  label 'biokanga'
  label 'align'
  tag {alignmeta}

  input:
    set val(simmeta), file("?.fq.gz"), val(dbmeta), file(kangadb) from reads4kangaAlign.combine(kangadbs) //cartesian product i.e. all input sets of reads vs all dbs

  output:
    set val(alignmeta), file('out.bam') into kangaBAMs

  when: //only align reads to the corresponding genome!
    simmeta.species == dbmeta.species && simmeta.version == dbmeta.version

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


// echo true
// categories = ["M_1":"1-st segment is correctly mapped", "M_2":"2-nd segment is correctly mapped",
//               "m":"segment should be unmapped but it is mapped", "w":"segment is mapped to a wrong location",
//               "U":"segment is unmapped and should be unmapped", "u":"segment is unmapped and should be mapped"]
// entries = Channel.create()
// summaries.subscribe { meta, f ->
//   entry = meta.clone()
//   entry.results = [:]
//   f.eachLine {  line ->
//     (k, v) = line.split()
//     entry.results << [(k) : v ]
//   }
//   //println entry.results
//   entries << entry
// }
// entries.subscribe{
//   println it
// }
// println(prettyPrint(toJson(entries)))




process collateResults {
  // echo true

  input:
    val collected from summaries.collect()
  
  output:
    file '*'

  exec:
  def outfile = task.workDir.resolve('results.json')
  categories = ["M_1":"1-st segment is correctly mapped", "M_2":"2-nd segment is correctly mapped",
  "m":"segment should be unmapped but it is mapped", "w":"segment is mapped to a wrong location",
  "U":"segment is unmapped and should be unmapped", "u":"segment is unmapped and should be mapped"]
  entry = null
  entries = []
  i=0;
  collected.each {
    if(i++ %2 == 0) {
      if(entry != null) {
        entries << entry
      }
      entry = it.clone()
    } else {
      entry.results = [:]
      // println "Current: $it"
      it.eachLine { line ->
        (k, v) = line.split()
        entry.results << [(categories[(k)]) : v ]
      }
    }
  }
  // println("FINAL: "+entries)
  outfile << prettyPrint(toJson(entries))

  // """
  // echo
  // """

}



// /*
//   Generic method for merging read meta with db meta and aligner info
//  */
//  def getAlignMeta(meta, dbmeta) {
//   return (meta + dbmeta) //.aligner =
// }