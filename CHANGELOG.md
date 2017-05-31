# SIMdrom Changelog

## development

### simdrom-cli

* Adapt API to new simdrom-core
* Introduce a seed option
* -- name to set a sample-name (overwrites name "Sampled" in vcf file)
* -bGC or -mGC
* het/hom/hemi support

### simdrom-core

* Two bugfixes in spike-in
  * Mutation was not inserted if the variant was the last in a contig
  * Error when two or mutations where inserted at the end of a contig
* big restructuring (using Builders, optionals, etc.)
* Fixing bug in ClinVarFilter (wrong object was used to write down the omim IDs)
* Fixing bugs in VCFSampler and SpikeIn
* Fixing seed option bugs (Collection.shuffle was still random even if a seed was set)
* Limiting possible alleles to max. 2. Something like 2/3/4 should not be possible anymore!
* fix in hardy weinberg
* add possibility to use hom/het/hemi allele counts from exac etc.
* Adding sex to simdrom
  * user can define the sex of the simulated individual. It can be also be randomly assigned or it can be none (only autosomes)
  * Hardy weinberg sampling works dependend on the sex on X and Y
  * Same for ac hom/get/hemi
* generate genotype not correctly (consider all possible genotype, not only 2, seen at the side when using )
* adding gnomAD support using the genotype counts (GC) flag


### manual

* creating readthedocs manual!

## v0.2

* no changelog available
