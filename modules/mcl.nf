process fna_mcl_clust{
  //directives
  //module "mcl"
  container "sysbiojfgg/mcl:v0.1"
  publishDir "plasx_predicition/mcl_clusters", mode: 'copy'
  cpus 4
  time 5.h

  input:
    path(distances)
  output:
    path("*.clusters"), optional: true, emit: clusters
  script:
    """
    if [ -s $distances ]; then
    mcl $distances --abc -te 4 -I $params.inflation -o ${distances}.clusters;
    else
    touch ${distances}.clusters
    fi
    """
}


process fna_get_rep{
  //directives
  publishDir "plasx_prediction", mode: 'copy'
  //conda "pandas"
  //module "python3"
  container "biopython/biopython"
  cpus 1
  time 4.h

  input:
    path(fna)
    path(clust)
  output:
    path("representative_plasmids.fna"), emit: representatives
  script:
    """
    get_representatives.py  --plasmids ${fna}  \\
                --clusters  ${clust} 
    """
}



