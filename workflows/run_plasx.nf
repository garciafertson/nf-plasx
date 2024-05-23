/*
PLasX: A workflow for the identification of plasmids in metagenomic data
*/

//import modules
include {anvio_prodigal} from  "../modules/anvio"
include {anvio_cogpfam} from  "../modules/anvio"
include {plasx_search_fam} from  "../modules/plasx"
include {plasx_predict} from  "../modules/plasx"
include {get_fna_plasmids} from  "../modules/plasx"
include {fna_mashtriangle} from  "../modules/mash"
include {fna_mcl_clust} from  "../modules/mcl"
include {fna_get_rep} from  "../modules/mcl"
include {rgi_card} from  "../modules/arg"
include {deep_arg} from  "../modules/arg"


workflow PLASX {
  //Define input parameters
  //Create channel contigs from fasta files
  Channel
    .fromPath(params.input)
    .set { ch_contigs }
  
  //Run anvio, predict genes with prodigal
  anvio_prodigal(ch_contigs)
  anvio_contigdb = anvio_prodigal.out.contigsdb
  genecalls = anvio_prodigal.out.genecalls
  allcontigfna = anvio_prodigal.out.fna.collect()

  //Predict COGS and Pfam v32
  anvio_cogpfam(anvio_contigdb)
  cogspfams = anvio_cogpfam.out.cogspfams

  // Use PlasX to search for search for de nove gene families
  plasx_search_fam(genecalls)
  plasxfams=plasx_search_fam.out.fams

  //Combine channels anvio_cogs_and_pfams, plasx_search_fams, and anvio_gene_calls by key
  cog_plasx= cogspfams.combine(plasxfams, by: 0)
  // .map { id, cog, plasxfam -> [cog, plasxfam] }
  
  //Combine channels genecalls, cog_plasx, and anvio_gene_calls by key
  cog_plasx_genecalls=cog_plasx.combine(genecalls, by: 0)
  //.map { id, cog, plasxfam -> [cog, plasxfam] }

  // Use PlasX to predict plasmids
  // and select  contigs where score of predicted plasmid is above threshold 0.7
  plasx_predict(cog_plasx_genecalls)
  plasmidsscores = plasx_predict.out.scores.collect()

  //////////////////////////////////////////
  // Create database of predicted plasmids
  // Combine predicted contigs
  get_fna_plasmids(allcontigfna, plasmidsscores)
  plasmidsfna = filter_fna_plasmids.out.plasmidsfna

  // Dereplicate plasmids using mash, 
  // cluster using mcl 
  // get representative sequences (longest sequence in cluster)
  fna_mashtriangle(plasmidsfna)
  distances = fna_mashtriangle.out.list05
  fna_mcl_clust(distances)
  clusters=fna_mcl_clust.out.clusters  //publish clusters to draw network
  fna_get_rep(df,asdb_fnafiles, asdb_distances, clusters)
  plasmid_representatives = fna_get_repr.out.representatives

  //anotate predicted plasmids with deepARG and rgi card database
  //retrieve annotation as contig and list of detected CARD genes
  rgi_card(plasmid_representatives)
  deep_arg(plasmid_representatives)
  
  //Combine with Plasmid database human samples,
  // 1 Blastn against plasmid database
  // 2 Get best hit sequences, long enough (at leas 50% of query length)
  // 3 Extract fasta sequence and annotations of blast hits in plasmid database
  // 
  
  //////////////////////////////////////////
  //Create index for bowtie2 alingment generate bam files
  

  // map reads to plasmid
 
 // Use  MobMess for inferring and visualizing  evolutionary relations among plasmid sequences
 //infer non redundant subset
 //infer plasmid systems, a backbone plasmid with core genes acceres accesory genes to form compoud plasmids
 //visualize similarity network of many plasmids
 //visualize plasmid systems, shared get and backbone genes of small set of plasmids
}
