* [![Build Status](https://travis-ci.org/visze/simdrom.svg?branch=master)](https://travis-ci.org/visze/simdrom) - Master 
* [![Build Status](https://travis-ci.org/visze/simdrom.svg?branch=development)](https://travis-ci.org/visze/simdrom) - Development
* [![Documentation Status](https://readthedocs.org/projects/simdrom/badge/?version=latest)](http://simdrom.readthedocs.io/en/latest/?badge=latest)
* [![Documentation Status](https://readthedocs.org/projects/simdrom/badge/?version=development)](http://simdrom.readthedocs.io/en/development/?badge=development)

# SIMdrom

Simulate Exomes, Genomes, or self defined populations and spike in mutations. SIMdrom is a Java program to randomly select variants in a VCF file and create a new VCF file. A second VCF file can be used to spike in other randomly selected variants. Homozygous and heterozygous genotypes are selected by the Hardy-Weinberg principle. SIMdrom also exposes its functionality through a library API.

SIMdrom is compatible with Java 8 and higher.

## Use cases

There are two typical use-cases where you can use SIMdrom. 

1. Generate a random variant files to test you NGS pipeline.
2. Generate a random variant files and spike in a pathogenic mutations. You can use the files for performance validation of your diagnostic, gene prioritisation, or variant prioritisation pipeline.   
 

## Quickstart

### Use 1000 Genomes to generate a random genome

Download the *.sites.vcf.gz from 1000Genomes - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ - and the corresponding index file *.sites.vcf.gz.tbi. The VCF-file should contain a `AF` identifier in the info column.

Run SIMdrom with the downloaded VCF file
```
# java -jar simdrom-cli-0.0.1.jar -b ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF AF
```
The new sampled individual will be printed out into the standard output. You can pipe it into a VCF file.
 ```
# java -jar simdrom-cli-0.0.1.jar -b ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF AF > newIndividualVCFfile.vcf
```
or generate a bgziped file
```
# java -jar simdrom-cli-0.0.1.jar -b ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF AF | bgzip -c | > newIndividualVCFfile.vcf.gz
```
or use the `--output` flag to print it to a (bgzip) file directly.
```
# java -jar simdrom-cli-0.0.1.jar -b ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF AF --output newIndividualVCFfile.vcf.gz
```

### Use ExAC to generate a random Exome

Download the current ExAC release from the Broadinstitute FTP - ftp://ftp.broadinstitute.org/pub/ExAC_release/ - and the corresponding index file. The VCF-file should contain a `AF` identifier in the info column.

Run SIMdrom with the downloaded VCF file
```
# java -jar simdrom-cli-0.0.1.jar -b ExAC.r0.3.sites.vep.vcf.gz -bAF AF
```

### Use 1000 Genomes to generate a random genome of an ethnical population

Download the 1000Genomes sites VCF file as described earlier. Ít contains allele frequencies from 5 different populations:

|Identifier | Population|
|-----------|-----------|
|EAS_AF     |EAS        |
|EUR_AF     |EUR        |
|AFR_AF     |AFR        |
|AMR_AF     |AMR        |
|SAS_AF     |SAS        |

For example if you want to sample an european individual use the following command:
```
# java -jar simdrom-cli-0.0.1.jar -b ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF EUR_AF
```

### Use ExAC to generate a random Exome of an ethnical population

Download the ExAC VCF file as described earlier. It contains allele frequencies from 7 different populations. But instead of allele frequencies the counts of alternative and total alleles are provided. So we have to use two options in SIMdrom `-bAC` (alt allele counts) and `-bAN`(total allele counts)



|Identifier ALT allele counts | Identifier total allele counts | Population           |
|-----------------------------|--------------------------------|----------------------|
|AC_AFR                       |AN_AFR                       |African/African American |
|AC_AMR                       |AN_AMR                       |American                 |
|AC_EAS                       |AN_EAS                       |East Asian               |
|AC_FIN                       |AN_FIN                       |Finnish                  |
|AC_NFE                       |AN_NFE                       |Non-Finnish European     |
|AC_OTH                       |AN_OTH                       |Other Allele             |
|AC_SAS                       |AN_SAS                       |South Asian Allele       |

For example if you want to sample an american individual use the following command:
```
# java -jar simdrom-cli-0.0.1.jar -b ExAC.r0.3.sites.vep.vcf.gz -bAC AC_AMR -bAN AN_AMR
```

### Randomly select one genome of 1000Genomes

You can also just select one genotype of an individual of a 1000genomes sample. Therefore you have to download the genotype VCF files from 1000Genomes - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ - and the corresponding index files. Right now, SIMdrom can only use one input VCF as background population. Therefore you have to merge the files that are divided by chromosome. Now you can use the `--single-sample` option to select only one genotype.
 ```
# java -jar ALL.ALLChr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --single-sample
```

Attach name of one sample after `--single-sample` and only this specific sample will be selected:
 ```
# java -jar ALL.ALLChr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --single-sample HG00113
```

### Spike in a pathogenic mutation of ClinVar

Download the ClinVar VCF ind index file from the NCBI FTP - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/. Not every included variant is known as pathogenic. To use only the known pathogenic variants we have to use the info column filter of SIMdrom on the mutations file: `--mutations-info-filter`. The pathogenicity is decoded in the `CLNSIG` identifier with the number `5`. The corresponding SIMdrom option is `--mutations-info-filter CLNSIG=5`.

If we only want to spike in 1 mutation we have to use the `--mutations-variants-amount` with the value `1`. We can increase the value if we want more spiked in mutations. Every mutation has an equal probability to be chosen.

To find out the spiked in mutation(s), we can deliver a spike in log-file with the option `--spike-in-log`. The spike in log format is TSV and includes every spiked in mutation (TSV for better readability). 

In this example we use ExAC as background population and spike in one ClinVar mutation:
```
# java -jar simdrom-cli-0.0.1.jar -b ExAC.r0.3.sites.vep.vcf.gz -bAF AF -m clinvar.vcf.gz --mutations-info-filter CLINSIG=5 --mutations-variants-amount 1 --spike-in-log clinVarSpikeInLog.tsv
```
