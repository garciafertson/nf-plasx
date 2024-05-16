// main workflow assembly SRA using megahit
nextflow.enable.dsl=2
include {PLASX} from './workflows/run_plasx'
//run assembly pipeline
workflow NF_PLASX {
    PLASX()
}
//WORKFLOW: Execute a single named workflow for the pipeline
workflow {
    NF_PLASX ()
}
