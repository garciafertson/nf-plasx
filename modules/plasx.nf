
process plasx_search_fam{
  //set directives to control
  //scratch true
  cpus '6'
  memory '16 GB'
  time '8h'
  maxForks 40
  container "sysbiojfgg/plasx:v0.1" 
  errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  
  input:
    tuple val(x), path(genecalls)
  output:
    tuple val(x), path("${x}-de-novo-families.txt"), emit: fams
  script:

  """
   plasx search_de_novo_families \\
    -g ${x}-gene-calls.txt \\
    -o ${x}-de-novo-families.txt \\
    --threads 6 \\
    --splits 32 \\
    --overwrite
  """
  }

process plasx_predict{ 
    //set directives to control
    //scratch true
    cpus '2'
    memory '16 GB'
    time '4h'
    maxForks 60
    container "sysbiojfgg/plasx:v0.1" 
    errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
    maxRetries 3
    publishDir "${params.outdir}/mge/plasx/scores", mode: 'copy'
    
    input:
        tuple val(x), path(cogs), path(fams), path(genecalls)
    output:
        tuple val(x), path ("${x}-plasmid.txt") , emit: scores
    script:
    
    """
    plasx predict \\
    -a ${cogs} ${fams} \\
    -g ${genecalls} \\
    -o ${x}-scores.txt \\
    --overwrite

    filter_plasmid_scores.py \\
    --plasmids ${x}-scores.txt \\
    --threshold ${params.plamsmid_threshold} \\
    --output ${x}-plasmid.txt
    """
}

