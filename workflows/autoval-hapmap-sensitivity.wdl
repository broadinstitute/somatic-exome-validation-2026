### This is the sensitivity part of agnostic autovalidation

### The original wdl (including the comments) is written by David Benjamin <davidben@broadinstitute.org>,
### and modified by Junko Tsuji <jtsuji@broadinstitute.org> to accept MAF or VCF inputs from any callers.

#Conceptual Overview
# To measure sensitivity, we sequence a pool 5, 10 or 20 normal samples. The pool contains a variety of allele fractions,
# like a real tumor sample. These normal samples were also sequenced individually as part of HapMap, so we have "truth" vcfs
# for them.  We calculate sensitivity by comparing the calls to the truth data.  For a variety of reasons we don't call
# against a matched normal.


#Workflow Steps
# 1.  Restrict a huge HapMap wgs VCF to given lists of samples and intervals.
# 2.  Annotate this VCF using a bam sequenced from a pool derived from the given samples
# 3.  Convert MAF to VCF if the variant input is in MAF
# 4.  Compare calls to the truth data and output a table of true positives and false negatives along with
#     annotations from the truth VCF prepared in steps 1 and 2.

# Here we implement these steps, except for the subsampling in step 1, for several replicate bams each of 5-plex, 10-plex, and 20-plex mixtures
# e.g. 4 different 10-plex bams, 3 10-plex bams, and 7 20-plex bams.

import "https://api.firecloud.org/ga4gh/v1/tools/GP-TAG:hapmap-single-plex/versions/22/plain-WDL/descriptor" as single_plex

workflow HapmapSensitivityAllPlexes {
    File? intervals
    Boolean is_maf_input

    Int max_depth
    Array[Int] depth_bins
    Int depth_bin_width

    # bam and maf/vcf list for 5, 10, and 20 plex    
    Array[File] replicates_lists
    Array[String] prefixes
    String outbasename
    
    # truth vcfs for 5, 10, and 20 plex
    Array[File] preprocessed_truth 
    Array[File] preprocessed_truth_idx

    File python_script
    File r_script
    
    # python script for extracting major entries in MAF to avoid errors
    File petite_maf_python

    # runtime
    File? gatk_override
    String gatk_docker
    String python_docker
    String pdftk_docker

    # oncotator parameters
    File? onco_ds_tar_gz
    String? oncotator_docker

    Int disk_padding = 20

    scatter (n in range(length(replicates_lists))) {
        call single_plex.HapmapSensitivity {
              input:
                  gatk_override = gatk_override,
                  gatk_docker = gatk_docker,
                  python_docker = python_docker,
                  oncotator_docker = oncotator_docker,
                  onco_ds_tar_gz = onco_ds_tar_gz,
                  is_maf_input = is_maf_input,
                  intervals = intervals,
                  max_depth = max_depth,
                  depth_bins = depth_bins,
                  depth_bin_width = depth_bin_width,
                  replicate_list = replicates_lists[n],
                  preprocessed_hapmap = preprocessed_truth[n],
                  preprocessed_hapmap_idx = preprocessed_truth_idx[n],
                  prefix = prefixes[n],
                  python_script = python_script,
                  r_script = r_script,
                  petite_maf_python = petite_maf_python
        }
    }

    call single_plex.CombineTables as AllPlexSensitivityTable {
        input:
            gatk_docker = gatk_docker,
            input_tables = HapmapSensitivity.raw_sensitivity_table,
            prefix = "all_plex_sensitivity",
            disk_space = disk_padding
    }

    call single_plex.CombineTables as AllPlexPrecisionTable {
        input:
            gatk_docker = gatk_docker,
            input_tables = HapmapSensitivity.raw_precision_table,
            prefix = "all_plex_precision",
            disk_space = disk_padding
    }

    call single_plex.AnalyzeSensitivity as AllPlexSensitivity {
        input:
            python_docker = python_docker,
            input_table = AllPlexSensitivityTable.table,
            python_script = python_script,
            prefix = "all_plex",
            depth_bins = depth_bins,
            depth_bin_width = depth_bin_width,
            disk_space = ceil(size(AllPlexSensitivityTable.table, "GB")) + disk_padding
    }

    call single_plex.AnalyzePerformance as AllPlexPerformance {
        input:
            precision_table = AllPlexPrecisionTable.table,
            sensitivity_table = AllPlexSensitivityTable.table,
            r_script = r_script,
            prefix = "all_plex",
            depth_bins = depth_bins,
            depth_bin_width = depth_bin_width,
            max_depth = max_depth,
            disk_space = ceil(size(AllPlexPrecisionTable.table, "GB") + size(AllPlexSensitivityTable.table, "GB")) + disk_padding
    }
  
    call MergePlots as MergeSensitivityPlots { 
        input: 
            combined_plot1 = AllPlexSensitivity.snp_plot, 
            combined_plot2 = AllPlexSensitivity.indel_plot,
            plots1 = HapmapSensitivity.snp_plot,
            plots2 = HapmapSensitivity.indel_plot,
            plots3 = [],
            output_name = outbasename + ".sensitivity_plots.pdf",
            pdftk_docker = pdftk_docker
    }


    call MergePlots as MergePerformancePlots { 
        input: 
            combined_plot1 = AllPlexPerformance.expected_estimated_af_plot,
            combined_plot2 = AllPlexPerformance.fp_estimated_af_plot,
            combined_plot3 = AllPlexPerformance.depth_roc_plot,
            plots1 = HapmapSensitivity.expected_estimated_af_plot,
            plots2 = HapmapSensitivity.fp_estimated_af_plot,
            plots3 = HapmapSensitivity.depth_roc_plot,
            output_name = outbasename + ".performance_plots.pdf",
            pdftk_docker = pdftk_docker
    }
  
    call Concatenate as ConcatenateSummaries {
        input:
            inputs1 = HapmapSensitivity.summary,
            inputs2 = [],
            inputs3 = [],
            docker = python_docker,
            output_name = outbasename + ".summary.txt"
    }

    call Concatenate as ConcatenateTruthStats {
        input:
            inputs1 = HapmapSensitivity.truth_stats,
            inputs2 = [],
            inputs3 = [],
            docker = python_docker,
            output_name = outbasename + ".truth_stats.txt"
    }
  
    call Concatenate as ConcatenateRaw {
        input:
            inputs1 = [AllPlexSensitivity.snp_table, AllPlexSensitivity.indel_table],
            inputs2 = HapmapSensitivity.snp_table,
            inputs3 = HapmapSensitivity.indel_table,
            docker = python_docker,
            output_name = outbasename + ".raw_tables.txt"
    }
  
    call Concatenate as ConcatenateJaccard {
        input:
            inputs1 = HapmapSensitivity.snp_jaccard_table,
            inputs2 = HapmapSensitivity.indel_jaccard_table,
            inputs3 = [],
            docker = python_docker,
            output_name = outbasename + ".jaccard.txt"
    }
  
    output {
        File raw_tables = ConcatenateRaw.concatenated
        File summary = ConcatenateSummaries.concatenated 
        File sensitivity_plots = MergeSensitivityPlots.plots
        File performance_plots = MergePerformancePlots.plots
        File jaccard = ConcatenateJaccard.concatenated
        File truth_stats = ConcatenateTruthStats.concatenated

        Array[Array[File]] filter_analysis = HapmapSensitivity.filter_analysis
        Array[Array[File]] tpfn = HapmapSensitivity.tpfn
        Array[Array[File]] tpfn_idx = HapmapSensitivity.tpfn_idx
        Array[Array[File]] tpfp = HapmapSensitivity.tpfp
        Array[Array[File]] tpfp_idx = HapmapSensitivity.tpfp_idx
        Array[Array[File]] ftnfn = HapmapSensitivity.ftnfn
        Array[Array[File]] ftnfn_idx = HapmapSensitivity.ftnfn_idx
    }
}

task Concatenate {
    Array[File] inputs1
    Array[File] inputs2
    Array[File] inputs3
    String docker
    String output_name
    Int? disk_space

    command {
        touch result
        for file in ${sep=' ' inputs1} ${sep=' ' inputs2} ${sep=' ' inputs3}; do
            name=`basename $file`
            echo $name >> result
            cat $file >> result
            echo "" >> result
        done

        mv result ${output_name}
    }

    runtime {
        docker: docker
        disks: "local-disk " + select_first([disk_space, 50]) + " HDD"
        memory: "1 GB"
        preemptible: 2
    }
    
    output {
        File concatenated = "${output_name}"
    }

}

task MergePlots {
    File combined_plot1
    File combined_plot2
    File? combined_plot3
    Array[File] plots1
    Array[File] plots2
    Array[File] plots3
    String pdftk_docker
    String output_name
    Int? disk_space

    command {
        pdftk ${combined_plot1} ${combined_plot2} ${combined_plot3} ${sep=' ' plots1} ${sep=' ' plots2} ${sep=' ' plots3} cat output ${output_name}
    }

    runtime {
        docker: pdftk_docker
        disks: "local-disk " + select_first([disk_space, 50]) + " HDD"
        memory: "2 GB"
        preemptible: 2
    }
    
    output {
        File plots = "${output_name}"
    }
}
