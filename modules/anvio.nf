process anvio_prodigal {
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
    tuple val(x), path("${x.id}.db"), emit: contigsdb
    tuple val(x), path("${x.id}-gene-calls.txt"), emit: genecalls
  script:

  """
  # - The `-L 0` parameter ensures that contigs remain intact and aren't split
    anvi-gen-contigs-database -L 0 -T 6 --project-name ${x.id} -f ${contigs} -o ${x.id}.db

  # Export gene calls (including amino acid sequences) to text file
  anvi-export-gene-calls --gene-caller prodigal -c ${x.id}.db -o ${x.id}-gene-calls.txt

  """
}

process anvio_cogpfam {
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
    tuple val(x), path("${x.id}-cogs-and-pfams.txt "), emit: cogspfams
  script:

  """
  anvi-run-pfams -T 6 --pfam-data-dir Pfam_v32 -c ${x.id}.db
  #Export functions to text file
  anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c ${x.id}.db -o ${x.id}-cogs-and-pfams.txt
  """
} 