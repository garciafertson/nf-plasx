process cdhit {
  cpus '4'
  memory '8 GB'
  time '6h'
  maxForks 40
  container "nanozoo/cdhit:4.8.1--c697693"
  //errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
  publishDir "arg/catalogue", mode: 'copy'
  
  input:
    val(name)
    tuple val(x), path(orfs)
  output:
    tuple val(x), path("${x}-cdhit.fna"), emit: gene_catalogue
    tuple val(x), path("${x}-cdhit.fna.clstr"), emit: clusters

  script:
    """
    cat ${orfs} > all_orfs.fna
    cd-hit -i all_orgs.fna -o ${name}-cdhit.fna -c 0.95 -n 5 -d 0 -M 8000 -T 4
    """
} 
