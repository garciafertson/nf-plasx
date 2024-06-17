process anvio_prodigal {
  cpus '4'
  memory '16 GB'
  time '4h'
  maxForks 40
  container "sysbiojfgg/anvio_cogpfam:v0.1" 
  errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  maxRetries 3
  
  input:
    path(contigs)
  output:
    tuple val(x), path("${x}.db"), emit: contigsdb
    tuple val(x), path("${x}-gene-calls.txt"), emit: genecalls
    tuple val(x), path("${x}-fixed.fa"), emit: fna
    tuple val(x), path("${x}-orfs.fna"), emit: orfs
  
  script:
  x=contigs.getSimpleName()
  """
  anvi-script-reformat-fasta ${contigs} -o ${x}-fixed.fa -l 1000 --simplify-names --prefix ${x}_
  # - The `-L 0` parameter ensures that contigs remain intact and aren't split
  anvi-gen-contigs-database -L 0 -T 4 --project-name ${x} -f ${x}-fixed.fa -o ${x}.db

  # Export gene calls (including amino acid sequences) to text file
  anvi-export-gene-calls --gene-caller prodigal -c ${x}.db -o ${x}-gene-calls.txt

  anvi-get-sequences-for-gene-calls -c ${x}.db -o ${x}-orfs.fna \\
  --report-extended-deflines

  #replace whitespace in fasta with hyphen
  sed -i 's/ /;/' ${x}-orfs.fna 
  
  """
}

process anvio_cogpfam {
  cpus '6'
  memory '16 GB'
  time '4h'
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

process anvio_hmm_mge {
cpus '6'
memory '8 GB'
time '4h'
input:
tuple val(x), path(contigsdb)
container "sysbiojfgg/anvio_cogpfam:v0.1"
containerOptions "--bind ${projectDir}/hmm_rec:${projectDir}/hmm_rec"

output:
tuple val(x), path("${x}_DOMTABLE.txt"), emit: hmm
tuple val(x), path("${x}_hmmrec.fa"), emit: rec_orfs

script:
"""
# Detect Recombinanses using HMM profiles
cp ${contigsdb} ${x}_hmm.db
anvi-run-hmms -c ${x}_hmm.db \\
              -H ${projectDir}/hmm_rec \\
              --num-threads 6 \\
              --hmmer-output-dir hmm-output \\
              --domain-hits-table

anvi-get-sequences-for-hmm-hits -c ${x}_hmm.db \\
                    --hmm-source hmm_rec \\
                    -o ${x}_hmmrec.fa

# mv domain hits to current directory
mv hmm-output/hmm.domtable ./${x}_DOMTABLE.txt
"""
}