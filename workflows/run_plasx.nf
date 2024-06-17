/*
PLasX: A workflow for the identification of plasmids in metagenomic data
*/
//import modules
include {anvio_prodigal} from  "../modules/anvio"
include {anvio_cogpfam} from  "../modules/anvio"
include {plasx_search_fam} from  "../modules/plasx"

include {plasx_predict} from  "../modules/plasx"
include {get_fna_plasmids} from  "../modules/biopython"

include {fna_mashtriangle as fna_mashtriangle_plasmid; fna_mashtriangle as fna_mashtriangle_rec} from  "../modules/mash"
include {fna_mcl_clust as fna_mcl_clust_plasmid; fna_mcl_clust as fna_mcl_clust_rec} from  "../modules/mcl"
include {fna_get_rep as fna_get_rep_plasmid; fna_get_rep as fna_get_rep_rec} from  "../modules/mcl"

include {rgi_card as rgi_card_sample} from  "../modules/arg"
include {deep_arg as deep_arg_sample; deep_arg as deep_arg_catalogue} from  "../modules/arg"
include {get_deeparg_fna} from  "../modules/biopython"

include {cdhit as cdhit; cdhit as cdhit_arg} from  "../modules/gene_catalogue"
include {anvio_hmm_mge} from  "../modules/anvio"
include {get_contig_from_hmm} from  "../modules/biopython"

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
  contigs = anvio_prodigal.out.fna
  orfs = anvio_prodigal.out.orfs
  
  //////////////////////////////////////////
  //RUN MGE PREDICTION
  if (params.run_mge_prediction) {
    ///////////////////////////////////
    ///////PLASMID PREDICTION
    ///////////////////////////////////
    // Use PlasX to search for search for de nove gene families
    plasx_search_fam(genecalls)
    plasxfams=plasx_search_fam.out.fams
    //Predict COGS and Pfam v32
    anvio_cogpfam(anvio_contigdb)
    cogspfams = anvio_cogpfam.out.cogspfams
    //Combine channels anvio_cogs_and_pfams, plasx_search_fams, and anvio_gene_calls by key
    cog_plasx= cogspfams.combine(plasxfams, by: 0)
    // .map { id, cog, plasxfam -> [cog, plasxfam] }
    //Combine channels genecalls, cog_plasx, and anvio_gene_calls by key
    cog_plasx_genecalls=cog_plasx.combine(genecalls, by: 0)
    // Use PlasX to predict plasmids
    // and select  contigs where score of predicted plasmid is above threshold 0.7
    plasx_predict(cog_plasx_genecalls)
    //Return annotation table with contig name, plasmid score
    plasmidsscores = plasx_predict.out.scores
    //// Create CATALOGUE of predicted plasmids
    // Combine predicted contigs
    contig_plasmidscore=contigs.combine(plasmidsscores, by: 0)
    get_fna_plasmids(contig_plasmidscore)
    plasmidsfna = get_fna_plasmids.out.plasmidsfna.collect()
    plasmid_contigs= tuple("PlasXcontigs", plasmidsfna)
    // Dereplicate plasmids using mash, 
    fna_mashtriangle_plasmid(plasmidsfna)
    distances = fna_mashtriangle_plasmid.out.list05
    allplasmids = fna_mashtriangle_plasmid.out.allplasmids
    // cluster using mcl 
    fna_mcl_clust_plasmid(distances)
    clusters=fna_mcl_clust_plasmid.out.clusters  //publish clusters to draw network
    // get representative sequences (longest sequence in cluster)
    fna_get_rep_plasmid(allplasmids, clusters)
    plasmid_representatives = fna_get_rep_plasmid.out.representatives 

    //////////////////////////////////////////
    //MGE RECOMBINASE PREDICTION
    //////////////////////////////////////////
    //Run hmmsearch against recombinase database
    anvio_hmm_mge(anvio_contigdb)
    //Returns fna orfs sequences
    mge_orfs= anvio_hmm_mge.out.rec_orfs
    mgeorg_contigs= mge_orfs.combine(contigs, by: 0)
    //Build orf hmm_mge table and recover contig names in text file
    //Subset contigs with recombinase hits
    get_contig_from_hmm(mgeorf_contigs)
    mge_contig = get_contig_from_hmm.out.contig_mge.collect()
    mge_contig = tuple("MGEcontigs", mge_contig)

    //Get representatives contigs with MGE hits
    fna_mashtriangle_rec(mge_contig)
    distances_reccontigs = fna_mashtriangle_rec.out.list05
    allplasmids_reccontigs = fna_mashtriangle_rec.out.allplasmids
    fna_mcl_clust_rec(distances_reccontigs)
    clusters_reccontigs=fna_mcl_clust_rec.out.clusters  
    fna_get_rep_rec(allplasmids_reccontigs, clusters_reccontigs)
    reccontigs_representatives = fna_get_rep_rec.out.representatives
  }
  //anotate predicted contigs with deepARG and rgi card database
  //retrieve annotation as contig and list of detected CARD genes
  //ANNOTATE CONTIGS WHERE ARGs ARE DETECTED
  if (params.run_ARG_prediction) {
    rgi_card_sample(contigs)
    // **** get rgi fasta sequence of SNIP hits and build ARG annotation table
    /*rgi_result = rgi_card_sample.out.rgi_fna
    rgi_getSNPseqs(rgi_result)
    rgi_SNPfna= rgi_getSNPseqs.out.rgi_SNPfna.collect()
    */

    deep_arg_sample(orfs)
    deepARG = deep_arg_sample.out.deeparg
    arg_orf= deepARG.combine(orfs, by: 0)
    get_deeparg_fna(arg_orf)
    deeparg_fna=get_deeparg_fna.out.deepargfna.collect()
    //Combine deeparg and rgi results and build ARG catalogue
    //arg_fna=deeparg_fna.concat(rgi_SNPfna)
    arg_fna= tuple ("ARGcatalogue", deeparg_fna)
    cdhit(arg_fna)

  }
  
  //Return annotation Contig NAME, ARG detected, MGE detected in all samples
  //RUN GENE CATALOGUE
  if(params.run_gene_catalogue) {
    //Build gene catalogue from predicted genes (prodigal) in assembly
    allorfs=orfs.collect()
    cdhit(allorfs)
    genecatlg= cdhit.out.fna
    //***** match ARG annotations from deepARG and rgi results in samples to gene catalogue
    //***** or run ARG prediction on gene catalogue
  }
}
