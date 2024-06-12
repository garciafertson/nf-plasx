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
  anvi-gen-contigs-database -L 0 -T 6 --project-name ${x} -f ${x}-fixed.fa -o ${x}.db

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

process anvio_hmm_mge {
cpus '6'
memory '8 GB'
time '4h'
input:
tuple val(x), path(contigsdb)

output:
tuple val(x), path("DOMTABLE_${x}.txt"), emit: hmm
tuple val(x), path("${x}_REC_contigs.fa"), emit: contig_fna

script:
"""
# Detect Recombinanses using HMM profiles
anvi-run-hmms -c ${contigsdb} \
              -H ${params.hmm_db} \
              --num-threads 6 \
              --hmmer-output-dir hmm-output \
              --domain-hits-table

# Filter HMM hits
anvi-script-filter-hmm-hits-table -c ${contigs-db} \
                                  --hmm-source ${params.hmm_db} \
                                  --domain-hits-table hmm-output/DOMTABLE.txt \
                                  --target-coverage 0.85

# Get fasta sequences of contigs hits
anvi-export-contigs -c ${contigs-db} \
                    -o path/to/${x}_REC_contigs.fa

# mv domain hits to current directory
mv hmm-output/DOMTABLE.txt ./DOMTABLE_${x}.txt
"""
}