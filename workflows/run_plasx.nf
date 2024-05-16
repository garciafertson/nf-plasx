/*
PLasX: A workflow for the identification of plasmids in metagenomic data
*/

//import modules
include {anvio_prodigal} from  "../modules/anvio"
include {anvio_cogpfam} from  "../modules/anvio"
include {plasx_search_fam} from  "../modules/plasx"
include {plasx_predict} from  "../modules/plasx"

workflow PLASX {
  //Define input parameters
  //Create channel contigs from fasta files
  Channel
    .fromPath(params.input)
    .ifEmpty { exit 1, "Cannot find any files matching: ${params.input}\nNB: Path needs to be enclosed in quotes!" }
    .map { row -> 
                def meta= [:]
                meta.id= row[0]
                return[meta, row[1]] 
                }
    .set { ch_contigs }
  
  //Run anvio, predict genes with prodigal
  anvio_prodigal(ch_contigs)
  anvio_contigdb = anvio_prodigal.out.contigsdb
  genecalls = anvio_prodigal.out.genecalls

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
  plasx_predict(cog_plasx_genecalls)
  plasx_predict.out.scores

 // Use  MobMess for inferring and visualizing  evolutionary relations among plasmid sequences
 //infer non redundant subset
 //infer plasmid systems, a backbone plasmid with core genes acceres accesory genes to form compoud plasmids
 //visualize similarity network of many plasmids
 //visualize plasmid systems, shared get and backbone genes of small set of plasmids
}
