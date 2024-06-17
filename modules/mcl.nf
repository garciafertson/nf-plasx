process fna_mcl_clust{
  //directives
  //module "mcl"
  container "sysbiojfgg/mcl:v0.1"
  publishDir "mge/catalogue/mcl_clusters", mode: 'copy'
  cpus 4
  memory "4 GB"
  time 5.h

  input:
    tuple val(x), path(distances)
  output:
    tuple val(x), path("*.clusters"), emit: clusters
  script:
    """
    if [ -s $distances ]; then
    mcl $distances --abc -te 4 -I $params.inflation -o ${x}_${distances}.clusters;
    else
    touch ${x}_${distances}.clusters
    fi
    """
}


process fna_get_rep{
  //directives
  publishDir "mge/catalogue/", mode: 'copy'
  container "biopython/biopython"
  memory "4 GB"
  cpus 1
  time 4.h

  input:
    path(fna)
    tuple val(x), path(clust)
  output:
    path("${x}_reprContigs.fna"), emit: representatives
  script:
    """
    get_representatives.py  --plasmids ${fna}  \\
                --clusters  ${clust}
    mv representative_plasmids.fna ${x}_reprContigs.fna 
    """
}

