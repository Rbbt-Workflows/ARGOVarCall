# Normal installation

If you had Rbbt installed the command would be

`rbbt workflow task VariantConsensus consensus_vcf --vcf_list <caller1.vcf>,<caller2.vcf>,<caller3.vcf> --caller_names <caller1>,<caller2>,<caller3>`

Where the vcfs are separated by commas and the names of the callers are also specified in the same order, or it will use the basenames of the files and it could get messy. You can also see more parameters doing

`rbbt workflow task VariantConsensus consensus_vcf -h`

# Singularity
 
I think the simplest way to use it is with singularity. There is an image here https://b2drop.bsc.es/index.php/s/3xJ9jmrz9RkQ4fL called rbbt-VariantConsensus.sif

Download it and do:

`singularity run rbbt-VariantConsensus.sif  rbbt workflow task VariantConsensus consensus_vcf --vcf_list <caller1.vcf>,<caller2.vcf>,<caller3.vcf> --caller_names <caller1>,<caller2>,<caller3>`

and it should work seamlessly like a regular install

# Docker

There is a also a docker container here https://hub.docker.com/repository/docker/mikisvaz/rbbt-variantconsensus

So you would do

`docker run <mounts> mikisvaz/rbbt-argovarcall rbbt workflow task VariantConsensus consensus_vcf --vcf_list <caller1.vcf>,<caller2.vcf>,<caller3.vcf> --caller_names <caller1>,<caller2>,<caller3>`

The problem with docker mounting the directories. The <mounts> is to make sure that the files from the callers are available to docker. Also, the intermediate files, which sill be created under '/root/.rbbt/var/jobs'. You can also change the directory that rbbt creates the files with

`docker run <mounts> mikisvaz/rbbt-argovarcall rbbt workflow task VariantConsensus consensus_vcf --vcf_list <caller1.vcf>,<caller2.vcf>,<caller3.vcf> --caller_names <caller1>,<caller2>,<caller3> --workdir_all <somedir>`

# Nextflow

It might be simpler to use the nextflow script that I've attached, which if you know nextflow you know you can call it with docker. There it assumes that all the vcfs are in a directory called 'data' and it doesn't use the --caller_names parameter, so the files better be called data/mutect2.vcf.gz, data/lofreq.vcf.gz, etc. Nextflow at least will take care of mounting stuff for you. Also, I've made it so that the nextflow step serves does the job in two steps, first to produced the combined file with all the info from all callers and then the consensus file with the cleaned up version

Any way, we can see it in more detail in a couple of weeks.
