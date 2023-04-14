#!/bin/env nextflow
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.vcf_list = "data/*.vcf"

process combine_vcf{
    input:
        path vcf_files
    output: 
        path 'combined.vcf'
    """
rbbt workflow task ARGOVarCall combine_vcf --vcf_list "$vcf_files" -O combined.vcf
    """
}

process consensus_vcf{
    input:
        path 'combined.vcf'
    output: 
        path 'consensus.vcf'

    """
rbbt workflow task ARGOVarCall consensus_vcf  -O consensus.vcf --override_deps 'ARGOVarCall#combine_vcf=combined.vcf'
    """
}

workflow {
    myFileChannel = Channel.fromPath( params.vcf_list )
    consensus_vcf(combine_vcf(myFileChannel.toList()))
}
