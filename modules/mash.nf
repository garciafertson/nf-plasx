
process fna2fnamash{
  //directives
  publishDir "bgc_catalogue/tmp_fnamash"
  //module "bioinfo-tools: mash"
  container "staphb/mash:2.3"
  cpus 1
  time 30.m

  input:
    path(fna)
  output:
    path("*.msh"), emit: mash
    path("filename"), emit: filename
  script:
    simplename=fna[0].getSimpleName()
    name=fna[0].getName()
    //msh="bgcmash/"+name+".msh"
    """
    mash sketch -o $name $fna
    echo "./${name}.msh" > filename
    """
}


process faa2faamash{
  //directives
  publishDir "bgc_catalogue/tmp_faamash"
  //module "bioinfo-tools: mash"
  container "staphb/mash:2.3"
  cpus 1
  time 30.m

  input:
    path(faa)
    output:
    path("*.msh"), emit: mash
    path("aafile"), emit: filename
  script:
    simplename=faa[0].getSimpleName()
    name=faa[0].getName()
    """
    mash sketch -a -o $name $faa
    echo "$launchDir/bgc_catalogue/tmp_faamash/${name}.msh" > aafile
    """
}


process mashdistance{
  //directives
  publishDir "bgc_catalogue/tmp_mashdistnace"
  //module "bioinfo-tools: mash"
  container "staphb/mash"
  //maxforks=100
  time '1h'
  cpus 4

  input:
    tuple path(pred_mash), path(rep_mash_list)
    path(allbgcpred)
  output:
    path ("${name}.distance.txt"), optional:true,	emit: distance
    path ("notin_asdb") , optional:true, emit: notinasdb
    path ("in_asdb"), optional:true, emit: inasdb
  script:
    name=pred_mash.getBaseName()
    """
    mash dist -p $task.cpus -d $params.mashdistance $pred_mash $rep_mash_list \\
    > ${name}.distance.txt
    cut -f 1 ${name}.distance.txt | sort | uniq > in_asdb
    comm -23 $allbgcpred in_asdb > tmpfilename
    completename.py --prefix $launchDir/bgc_catalogue/tmp \\
      --filenames tmpfilename
    """
	}


process mashtriangle{
  //directives
  //module "bioinfo-tools:mash"
  container "staphb/mash"
  publishDir "bgc_catalogue/tmp_mashtriangle"
  cpus params.mashcores
  time 24.h

  input:
    path(mashlist)
    path(mashdir)
    output:
    path("*.edgelist"), emit: edgelist
    path("*_0.05.list"), emit: list05
  script:
      """
      mash triangle -p 10 -E  -d 0.3 -l $mashlist > ${mashlist}.edgelist
      awk '{if (\$3 < $params.mashdistance) print \$1,\$2,\$3}' \\
      ${mashlist}.edgelist > ${mashlist}_0.05.list
      """
}


process fnamash_list{
  //directives
  //module "bioinfo-tools:mash"
  container "staphb/mash"
  time 1.h
  cpus 4

  input:
    path(list)
  output:
    path("*msh") , emit: mash
  script:
    """
    mash sketch -p $task.cpus -l $list 2> /dev/null
    """
}

process faamash_list{
  //directivers
  //module "bioinfo-tools:mash"
  container "staphb/mash"
  time 1.h
  cpus 4

  input:
    path(list)
    output:
    path("*msh"), emit:mash
  script:
    """
    mash sketch -p $task.cpus -a -l $list 2> /dev/null
    """
}
