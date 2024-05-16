process anvio_prodigal {
  cpus '6'
  memory '16 GB'
  time '8h'
  maxForks 40
  container "sysbiojfgg/anvio_cogpfam:v0.1" 
  errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  
  input:
    path(contigs)
  output:
    tuple val(x), path("${x}.db"), emit: contigsdb
    tuple val(x), path("${x}-gene-calls.txt"), emit: genecalls

  script:
  x=contigs.getSimpleName()
  """
  anvi-script-reformat-fasta ${contigs} -o ${x}-fixed.fa -l 1000 --simplify-names
  # - The `-L 0` parameter ensures that contigs remain intact and aren't split
  anvi-gen-contigs-database -L 0 -T 6 --project-name ${x} -f ${x}-fixed.fa -o ${x}.db
  # Export gene calls (including amino acid sequences) to text file
  anvi-export-gene-calls --gene-caller prodigal -c ${x}.db -o ${x}-gene-calls.txt
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
    tuple val(x), path("${x}-cogs-and-pfams.txt"), emit: cogspfams
  script:

  """
  anvi-run-ncbi-cogs -T 6 --cog-version COG14 --cog-data-dir /home/COG_2014 -c ${contigs}
  anvi-run-pfams -T 6 --pfam-data-dir /home/Pfam_v32 -c ${contigs}
  #Export functions to text file
  anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c ${contigs} -o ${x}-cogs-and-pfams.txt
  """
} 