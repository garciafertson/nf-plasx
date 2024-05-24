
process fna_mashtriangle{
  //directives
  container "staphb/mash:2.3"
  cpus 10
  time 24.h
  publishDir "plasx_prediction/mash_distances", mode: "copy"

  input:
    path(fna)
  output:
    path("plasmids.edgelist"), emit: edgelist
    path("plasmid_0.05.list"), emit: list05
  script:
    """
    mash triangle -p 10 \\
    -i -E ${fna} -k 19 > plasmids.edgelist

    awk '{if (\$3 < $params.mashdistance) print \$1,\$2,\$3}' plasmids.edgelist > plasmid_0.05.list
    """
}
