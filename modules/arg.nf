process rgi_card {
  cpus '4'
  memory '16 GB'
  time '10h'
  maxForks 4
  container "quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0" 
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  publishDir "mge_arg", mode: 'copy'
  
  input:
    path(plasmids)
  output:
    path("rgi_predicted_ARG.txt")

  script
  """
  rgi main --input_sequence ${plasmids} \\
           --output_file predicted_ARG.txt \\
           --input_type contig --clean \\
           --num_threads 4 --data plasmid \\
           --low_quality \\
           --split_prodigal_jobs
  """
}

process deep_arg{
    cpus '4'
    memory '16 GB'
    time '10h'
    maxForks 4
    container "gaarangoa/deeparg:latest"
    //container "quay.io/biocontainers/deeparg:1.0--py_0" 
    maxRetries 3
    publishDir "mge_arg", mode: 'copy'
    
    input:
        path(plasmids)
    output:
        path("deepARG_output.ARG")
    
    script
    """
    deeparg predict \\
    --model LS \\
    -i ${plasmids} \\
    -o deepARG_output \\
    --type nucl \\
    --min-prob 0.8 \\
    --arg-alignment-identity 30 \\
    --arg-alignment-evalue 1e-10 \\
    --arg-num-alignments-per-entry 1000 \\
    # -d /root \\
 
    """
}