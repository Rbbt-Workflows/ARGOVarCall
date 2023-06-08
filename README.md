Calculate consensus somatic variant calls from a set of VCFs

This workflow takes a number of steps to build the consensus calls. The process
starts by generating a 'combined' VCF file that contains all the calls from all
the callers. When two callers have called the same variant their information is
joined. To avoid overlaps in the INFO and FORMAT fields all the original fields
names are prefixed with the name of the caller (e.g. `mutect2--`). These new
names are added to the VCF header to keep the format valid. Since many tools
expect several fields to be present, like AD, AF, DP, etc., a short selection
of basic fields is extracted without prefixing, the actual values for these
fields are taken from the first caller that included them, in the order they
were specified as inputs. Note that this 'combined' VCF has contains all the
information from the original VCF files, specifically all the INFO and FORMAT
fields, but also including the file preambles (the commented lines at the
start), after removing repetitions.

The rest of the process finds out InDel variants that are called by the callers
in different ways. This is done by identifying sets of mutations that affect
overlapping segments of the sequence and then resolving what each caller thinks
the resulting sequence would be. When two callers agree on the resulting
sequence, the corresponding mutations are considered as 'equivalent' (when
considered collectively). To avoid counting the same mutational events multiple
times, each set of equivalent mutations sequence is assigned a
'UniqueRegionID'. The set of equivalent mutations is used to annotated the
'combined' VCF with the fields: CalledBy, ValidatedBy and UniqueRegionID. The
field CalledBy lists the callers that have called each individual mutation
exactly like that, while the ValidatedBy also include callers that have called
mutations in the list of mutations equivalent to the current one. Two
additional fields count the number of callers in each list: NumCallers,
NumValidators. 

Once all the information in the 'combined' VCF is ready the last step is to
fill in the FILTER field with FAIL or PASS depending on if number of callers in
'ValidatedBy' (or optionally CalledBy) is larger than the specified number. To
make the VCF more practical for uses a cleaned up version of it is produced at
the end, where all the prefixed fields are removed, leaving only the basic
common fields described earlier.

# Tasks

## combine_vcf
VCF file containing all the information from the individual caller VCF files

The information from each caller is including prefixing the name of the field with the name
of the caller. The basic common FORMAT fields (e.g. GT, AF, DP, AU, and AD) are taken
from the first caller that reports them.

## mutation_positions           

Extract the position of each mutation

## ranges

Turn the positions into ranges, by considering the length of InDel changes

## overlaps

Find regions of sequence that have sets of overlapping mutations

## reconstructed_sequence

For each overlap region, calculate the resulting mutated sequence for each caller

## equivalent_changes

For each overlap region, find callers that agree on the mutated sequence

## equivalent_mutations          

Form equivalence sets with mutations from callers that agree on each region

For each mutation it lists the callers that have called it, the UniqueRegionID
that contains it, and other mutations that collectively lead to the same
consequence in that UniqueRegionID. The lists of callers that concur on that
consequence is also listed.

## extended_validation_combined_vcf

Extend the combined VCF with the information about equivalent mutations

When a mutation is part of an equivalence set that includes a caller, this
caller is considered to also validate that mutation.

## fully_annotated_consensus_vcf

Consensus VCF where the FILTER field is field considering the number of callers that validate each mutation

## consensus_vcf

Remove original caller fields from consensus VCF

