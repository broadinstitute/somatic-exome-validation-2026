## Copyright Broad Institute, 2017
##
## This WDL workflow counts the number of variants (i.e. false positives) and
## reports false positive rates.
##
## Main requirements/expectations :
## - False Positive VCF or MAF files from normal vs normal calls
##
## Outputs :
## - Summary of false positive calls
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## pages at https://hub.docker.com/r/broadinstitute/* for detailed licensing information
## pertaining to the included programs.

## Modified by Junko Tsuji <jtsuji@broadinstitute.org>


import "https://api.firecloud.org/ga4gh/v1/tools/GP-TAG:maf-to-vcf/versions/10/plain-WDL/descriptor" as maf_convert

workflow SpecificityNormalNormal {
    Boolean is_maf_input
    File? intervals
    Array[File] variants
    Array[File]? variants_index
    String outbasename

    File? gatk_override
    String gatk_docker
    Int? preemptible_attempts

    File ref_fasta
    File ref_fai
    File ref_dict
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fai, "GB") + size(ref_dict, "GB"))
    Int disk_padding = 20

    # oncotator parameters
    File? onco_ds_tar_gz
    String? oncotator_docker
    
    # python script for extracting major entries in MAF to avoid errors
    File? petite_maf_python

    if (is_maf_input) {
        scatter (variant in variants) {
            call maf_convert.MafToVcfConversion {
                input:
                    input_maf = variant,
                    petite_maf_python = petite_maf_python,
                    sample_id = basename(variant, ".maf"),
                    onco_ds_tar_gz = onco_ds_tar_gz,
                    oncotator_docker = oncotator_docker,
                    gatk_docker = gatk_docker,
                    gatk_override = gatk_override
            }
        }
    }
    
    Array[File] vcf = select_first([MafToVcfConversion.output_vcf, variants])
    Array[File] vcf_index = select_first([MafToVcfConversion.output_vcf_idx, variants_index])
    scatter (n in range(length(vcf))) {
        call CountFalsePositives {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                input_vcf = vcf[n],
                input_vcf_index = vcf_index[n],
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = ref_size + ceil(size(vcf[n], "GB") + size(vcf_index[n], "GB")) + disk_padding
        }
    }

    call GatherTables {
        input:
            tables = CountFalsePositives.false_positive_counts,
            gatk_docker = gatk_docker,
            output_name = outbasename + ".summary.txt"
    }

    output {
        File summary = GatherTables.summary
    }
}

task GatherTables {
    # we assume that each table consists of two lines: one header line and one record
    Array[File] tables
    String gatk_docker
    String output_name

    command {
        set -e

        # extract the header from one of the files
        head -n 1 `ls ${sep=" " tables} | tr " " "\n" | head -n 1` > ${output_name}

        # then append the record from each table
        for table in `echo ${sep=" " tables}`; do
            tail -n +2 $table >> ${output_name}
        done
    }

    runtime {
        docker: gatk_docker
        memory: "1 GB"
        disks: "local-disk 20 HDD"
    }
    output {
        File summary = "${output_name}"
    }
}

task CountFalsePositives {
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_vcf
    File input_vcf_index

    File? gatk_override
    String gatk_docker

    Int? disk_space

    String vcf_basename = basename(input_vcf)

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        # move VCF and VCF index in the current directory
        mv ${input_vcf} ${input_vcf_index} .

        gatk --java-options "-Xmx4g" CountFalsePositives \
             -V ${vcf_basename} \
             -R ${ref_fasta} \
             ${"-L " + intervals} \
             -O "false-positives.txt"
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: "5 GB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: 3
    }

    output {
        File false_positive_counts = "false-positives.txt"
    }
}
