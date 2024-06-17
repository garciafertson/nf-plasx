process get_contig_from_hmm{
    //directives
    publishDir "mge/hmmprediction/table_seqs", mode: 'copy'
    container "biopython/biopython"
    cpus 1
    time 4.h

    input:
        tuple val(x), path(contigs), path(hmm) 
    output:
        path("${x}_recContig.fna"), emit: mge_contigs
        path("${x}_rechmm.fna"), emit: mge_orfs
        path("${x}_rechmm.tsv"), emit: mge_table

    script:
    """
    contig_from_anviohmm.py --hmm ${hmm} --output ${x} --contigs ${contigs}
    """
}

process get_fna_plasmids{
  cpus '1'
  memory '8 GB'
  time '1h'
  //maxForks 4
  container "biopython/biopython"
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 2
  publishDir "mge/plasx/plasmid_seqs", mode: 'copy'
  
  input:
    tuple val(x), path(contigfna) , path(plasmidscores)
  output:
    path("${x}_predicted_plasmids.fna"), emit: plasmidsfna
  script:
  """
  get_plasmid_fna.py --contigs ${contigfna} \\
  --plasmids ${plasmidscores} \\
  --output ${x}_predicted_plasmids.fna
  """
}