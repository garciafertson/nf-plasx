process gene_catalogue {
  cpus '6'
  memory '16 GB'
  time '8h'
  maxForks 40
  container "sysbiojfgg/anvio_cogpfam:v0.1" 
  errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  
  input:
    tuple val(x), path(contigs)
  output:
    tuple val(x), path("${x}-cogs-and-pfams.txt"), emit: cogspfams
  script:

  """
  anvi-run-ncbi-cogs -T 6 --cog-version COG14 --cog-data-dir /home/COG_2014 -c ${contigs}
  anvi-run-pfams -T 6 --pfam-data-dir /home/Pfam_v32 -c ${contigs}
  #Export functions to text file
  anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c ${contigs} -o ${x}-cogs-and-pfams.txt
  """
} 
