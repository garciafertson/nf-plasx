process mcl_clust{
  //directives
  //module "mcl"
  container "sysbiojfgg/mcl:v0.1"
  publishDir "bgc_catalogue/tmp_mashtriangle"
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


process get_rep_asdb{
  //directives
  publishDir "bgc_catalogue", mode: 'copy'
  //conda "pandas"
  //module "python3"
  container "biopython/biopython"
  cpus 1
  time 4.h

  input:
    path(df)
    path(mashlist)
    path(distedges)
    path(clust)
  output:
    path("*.tsv"), emit: tsv
    path("*.list"), emit: list
  script:
    """
    get_representatives.py  --asdbtsv $df \\
                --allfiles $mashlist \\
                --distance $distedges \\
                --clusters $clust \\
                --ext $params.ext \\
                --launchDir $launchDir
    """
}


process get_rep_predicted{
  //directives
  publishDir "bgc_catalogue", mode: 'copy'
  //conda "pandas"
  //module "python3"
  container "biopython/biopython"
  cpus 1
  time 1.h

  input:
    path(mashlist)
    path(distedges)
    path(clust)
    output:
    //path("predictedbgc_representative.tsv"), emit:tsv
    path("*.list"), emit:list

  script:
  """
  get_representatives_frompred.py --allfiles $mashlist \\
                                --distance $distedges \\
                                --clusters $clust \\
  """
}



process get_rep_predinasdb{
  //directives
  publishDir "bgc_catalogue", mode: 'copy'
  cpus 1
  //module "python3"
  container "biopython/biopython"
  time 1.h

  input:
    path(dist)
  output:
    path("*.list"), emit: list
    path("*.dict"), emit: dict
  script:
    """
    get_representatives_predinasdb.py --paired_distances $dist \\
                                    --ext $params.ext
    """
}
