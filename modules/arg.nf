process rgi_card {
  cpus '4'
  memory '16 GB'
  time '8h'
  maxForks 4
  container "quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0"
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  publishDir "arg/rgi", mode: 'copy'
  
  input:
    tuple val(x), path(contigs)
  output:
    tuple val(x), path("${x}_rgi_predicted_ARG.txt"), emit: rgi

  script:
  """
  rgi main --input_sequence ${contigs} \\
           --output_file ${x}_rgi_predicted_ARG \\
           --input_type contig --clean \\
           --num_threads 4 --data plasmid \\
           --low_quality \\
           --split_prodigal_jobs
  """
}

process deep_arg{
    cpus '2'
    memory '12 GB'
    time '8h'
    maxForks 4
    container "gaarangoa/deeparg:latest"
    containerOptions "--bind ${params.deeparg_db},${params.user_home}"
    //container "quay.io/biocontainers/deepargls:1.0--py_0" 
    maxRetries 3
    publishDir "arg/deeparg", mode: 'copy'
    
    input:
      tuple val(x),  path(orfs)
    output:
      tuple val(x),  path("${x}_deepARG.mapping.ARG"), emit: deeparg
    
    script:
    """
    deeparg predict \\
    --model LS \\
    -i ${orfs} \\
    -d ${params.deeparg_db} \\
    -o ${x}_deepARG \\
    --type nucl \\
    --min-prob 0.8 \\
    --arg-alignment-identity 30 \\
    --arg-alignment-evalue 1e-10 \\
    --arg-num-alignments-per-entry 1000 \\
   
    """
}