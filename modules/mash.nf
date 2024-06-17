
process fna_mashtriangle{
  //directives
  container "staphb/mash:2.3"
  cpus 10
  time 24.h
  memory "32 GB"
  publishDir "mge/catalogue/mashdistance", mode: "copy"

  input:
    val(x)
    path(fna)
  output:
    path("${x}_allplasmids.fna"), emit: allplasmids
    tuple val(x), path("${x}_plasmids.edgelist"), emit: edgelist
    tuple val(x), path("${x}_plasmid_0.05.list"), emit: list05
  script:
    """
    #concatenate plasmids
    cat *.fna > ${x}_allplasmids.fna

    mash triangle -p 10 \\
    -i -E ${x}_allplasmids.fna -k 19 > ${x}_plasmids.edgelist

    awk '{if (\$3 < $params.mashdistance) print \$1,\$2,\$3}' ${x}_plasmids.edgelist > ${x}_plasmid_0.05.list
    """
}
