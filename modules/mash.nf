
process fna_mashtriangle{
  //directives
  container "staphb/mash:2.3"
  cpus 10
  time 24.h
  memory "32 GB"
  publishDir "${params.outdir}/mge/catalogue/mashdistance", mode: "copy"

  input:
    val(x)
    path(fna)
  output:
    path("${x}_all.fna"), emit: allplasmids
    tuple val(x), path("${x}_all.edgelist"), emit: edgelist
    tuple val(x), path("${x}_0.05.list"), emit: list05
  script:
    """
    #concatenate plasmids
    cat *.fna > ${x}_all.fna

    mash triangle -p 10 \\
    -i -E ${x}_all.fna -k 19 > ${x}_all.edgelist

    awk '{if (\$3 < $params.mashdistance) print \$1,\$2,\$3}' ${x}_all.edgelist > ${x}_0.05.list
    """
}
