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
referencesRemoteFasta.mix(referencesLocal).into{ references4rnfSimReads; references4kangaIndex; references4bwaIndex; references4bowtie2Index } //references4kangaSimReads;

// simulators = ["ArtIllumina", "DwgSim", "MasonIllumina", "WgSim"] //also CuReSim if SE only
// reads4kangaAlign = Channel.create()
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
    set val(simmeta), file("*.fq.gz") into reads4bwaAlign, reads4bowtie2align //,reads4kangaAlign


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
    snakemake && gzip --fast *.fq   
    """
}

// process kangaIndex {
//   label 'biokanga'
//   tag{tag}
//   input:
//     set val(meta), file(ref) from references4kangaIndex

//   output:
//     set val(dbmeta), file(kangadb) into kangadbs

//   script:
//     tag=getTagFromMeta(meta, ' ')
//     dbmeta = ["species": meta.species, "version": meta.version]
//     """
//     biokanga index \
//     -i ${ref} \
//     -o kangadb \
//     --ref ${ref}
//     """
// }

process bwaIndex {
  label 'bwa'
  tag{tag}
  input:
    set val(meta), file(ref) from references4bwaIndex

  output:
    set val(dbmeta), file("${ref}*") into bwadbs

  script:
    tag=getTagFromMeta(meta, ' ')
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    bwa index -a bwtsw ${ref}
    """
}

process bowtie2Index {
  label 'bowtie2'
  tag{tag}
  input:
    set val(meta), file(ref) from references4bowtie2Index

  output:
    set val(dbmeta), file("${basename}*") into bowtie2dbs

  script:
    tag=getTagFromMeta(meta, ' ')
    basename=getTagFromMeta(meta)
    dbmeta = ["species": meta.species, "version": meta.version]
    """
    bowtie2-build ${ref} ${basename}
    """
}

// // kangadbsFeedback = Channel.create()
// process kangaAlign {
//   label 'biokanga'
//   tag {alignmeta}
//   input:
//     each metaAndReads from reads4kangaAlign //tuple is: set val(meta0),val(simmeta),file("?.fq.gz") 
//     // set val(simmeta),file("?.fq.gz") from reads4kangaAlign//.mix(kangaFASTQ)
//     set val(dbmeta),file(kangadb) from kangadbs//.mix(kangadbsFeedback)

//   //  output:
//   //   set val(dbmeta),file(kangadb) into kangadbsFeedback
//   //   set val(meta), val(tag), file("${tag}.bam") into kangaBAMs

//   // when:
//   //   meta0.species == dbmeta.species

//   script:
//     // meta = meta0.clone() //otherwise modifying orginal map, triggering re-runs with -resume
//     // meta.ref = dbmeta
//     // meta.aligner = "BioKanga"    
//     simmeta = metaAndReads[0]
//     basename = simmeta.tag+"_vs_"+dbmeta+".biokanga"
//     // tag = simmeta.toString()+" VS "+dbmeta 
//     alignmeta = dbmeta + simmeta
//     if(simmeta.mode == "SE") {
//       r1 = file(metaAndReads[1])
//       """
//       biokanga align \
//       -i ${r1} \
//       --sfx ${kangadb} \
//       --threads ${task.cpus} \
//       -o "${basename}.bam" 
//       """
//     } else if(simmeta.mode == "PE"){
//       r1 = file(metaAndReads[1][0])
//       r2 = file(metaAndReads[1][1])
//       """
//       biokanga align \
//       -i ${r1} \
//       -u ${r1} \
//       --sfx ${kangadb} \
//       --threads ${task.cpus} \
//       -o "${basename}.bam" \
//       --pemode 2 
//       """
//     }
// }

process bwaAlign {
  label 'bwa'
  tag {alignmeta}

  input:    
    each metaAndReads from reads4bwaAlign //tuple is: set val(meta0),val(simmeta),file("?.fq.gz") 
    set val(dbmeta), file('*') from bwadbs

  script:
    dbBasename=getTagFromMeta(dbmeta)
    simmeta = metaAndReads[0]
    alignmeta = dbmeta + simmeta

    // if(metaAndReads[2].getClass() instanceof ArrayList) {    
    
    if(simmeta.mode == 'SE') {
      r1 = file(metaAndReads[1])
      """
      bwa mem ${dbBasename}.fasta ${r1} > out.sam
      """
    } else {
      r1 = file(metaAndReads[1][0])
      r2 = file(metaAndReads[1][1])
      """
      bwa mem ${dbBasename}.fasta ${r1} ${r2} > out.sam
      """
    }
      
}



// bowtie2dbsFeedback = Channel.create()
process bowtie2align {
  label 'bowtie2'
  label 'samtools'

  tag {alignmeta}
  input:
    // set val(simmeta), file("?.gz") from reads4bowtie2align
    each metaAndReads from  reads4bowtie2align //tuple is: set val(meta0),val(simmeta),file("?.fq.gz") 
    set val(dbmeta), file("*") from bowtie2dbs//.mix(bowtie2dbsFeedback)
  
  // output:
  //   set val(dbmeta), file("${dbBasename}.*") into bowtie2dbsFeedback //re-using same db multiple times

  //TODO:dbBasename=getTagFromMeta(dbmeta)D when: condition to only align reads to the corresponding genome!

  script:
  dbBasename=getTagFromMeta(dbmeta)
  simmeta = metaAndReads[0]
  alignmeta = dbmeta + simmeta //.groupBy { it.id }.collect { it.value.collectEntries { it } }
  if(simmeta.mode == 'SE') {
    r1 = file(metaAndReads[1])
    """
    bowtie2 -x ${dbBasename} -U ${r1} -p ${task.cpus} | samtools view -bS > out.bam
    """
  } else {
    r1 = file(metaAndReads[1][0])
    r2 = file(metaAndReads[1][1])
    """
    bowtie2 -x ${dbBasename} -1 ${r1} -2 ${r2} -p ${task.cpus} | samtools view -bS > out.bam
    """
  }
}

/* 
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

// /* 
//   Generic method for merging read meta with db meta and aligner info
//  */
//  def getAlignMeta(meta, dbmeta) {
//   return (meta + dbmeta) //.aligner = 
// }