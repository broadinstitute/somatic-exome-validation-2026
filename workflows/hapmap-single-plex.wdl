### This is the sensitivity part of agnostic autovalidation

### The original wdl (including the comments) is written by David Benjamin <davidben@broadinstitute.org>,
### and modified by Junko Tsuji <jtsuji@broadinstitute.org> to accept MAF or VCF inputs from any callers,
### and to restrict GT and GQ values in the truth set.


#Conceptual Overview
# To measure sensitivity, we sequence a pool 5, 10 or 20 normal samples. The pool contains a variety of allele fractions,
# like a real tumor sample. These normal samples were also sequenced individually as part of HapMap, so we have "truth" vcfs
# for them.  We calculate sensitivity by comparing the calls to the truth data.  For a variety of reasons we don't call
# against a matched normal.


#Workflow Steps
# 1.  Restrict a huge HapMap wgs VCF to given lists of samples and intervals
# 2.  Filter variants with low confidence GT and GQ in the restricted HapMap VCFs
# 3.  Annotate this VCF using a bam sequenced from a pool derived from the given samples
# 4.  Convert MAF to VCF if the variant intput is in MAF
# 5.  Compare calls to the truth data and output a table of true positives and false negatives along with
#     annotations from the truth VCF prepared in steps 1 and 2.

# Here we implement these steps, except for the subsampling in step 1, for several replicate bams of a *single* plex,
# e.g. 4 different 10-plex bams.  Note that each replicate has its own truth vcf because although the truth variants are identical,
# the bam-derived annotations differ slightly.


import "https://api.firecloud.org/ga4gh/v1/tools/GP-TAG:maf-to-vcf/versions/10/plain-WDL/descriptor" as maf_convert

workflow HapmapSensitivity {
    File? intervals

    # list of bams and the corresponding mafs or vcfs
    Boolean is_maf_input
    File replicate_list
    Array[Array[String]] replicates = read_tsv(replicate_list)

    # prefix is for output table
    String prefix
    File python_script
    Int max_depth
    Array[Int] depth_bins
    Int depth_bin_width
    File preprocessed_hapmap
    File preprocessed_hapmap_idx
    Int disk_padding = 20
    
    # TAG scripts
    # python script for extracting major entries in MAF to avoid errors
    File petite_maf_python
    # python script to restrict GQ and missing GT
    File genotype_filter_python
    Float? genotype_rate
    Int? genotype_quality_cutoff
    Int? read_depth_cutoff
    # R script to generate more performance statistics
    File r_script

    File? gatk_override
    String gatk_docker
    String python_docker

    # oncotator parameters
    File? onco_ds_tar_gz
    String? oncotator_docker


    Int preprocessed_hapmap_size = ceil(size(preprocessed_hapmap, "GB") + size(preprocessed_hapmap_idx, "GB"))

    call RestrictIntervals {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            vcf = preprocessed_hapmap,
            vcf_idx = preprocessed_hapmap_idx,
            intervals = intervals,
            genotype_filter_python = genotype_filter_python,
            genotype_rate = genotype_rate,
            read_depth_cutoff = read_depth_cutoff,
            genotype_quality_cutoff = genotype_quality_cutoff,
            disk_space = preprocessed_hapmap_size + disk_padding
    }

    Int truth_vcf_size = ceil(size(RestrictIntervals.output_vcf, "GB") + size(RestrictIntervals.output_vcf_idx, "GB"))

    scatter (row in replicates) {
        File replicate_id = row[0]
        File bam = row[1]
        File bam_index = row[2]
        File variant = row[3]
        File variant_index = row[4]

        Int bam_size = ceil(size(bam, "GB") + size(bam_index, "GB"))

        call MakeTruth {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                vcf = RestrictIntervals.output_vcf,
                vcf_idx = RestrictIntervals.output_vcf_idx,
                bam = bam,
                bam_idx = bam_index,
                max_depth = max_depth,
                disk_space = bam_size + truth_vcf_size + disk_padding
        }

        # convert MAF to VCF if the input is MAF
        if (is_maf_input) {
            call maf_convert.MafToVcfConversion {
                input:
                    input_maf = variant,
                    sample_id = replicate_id,
                    petite_maf_python = petite_maf_python,
                    onco_ds_tar_gz = onco_ds_tar_gz,
                    oncotator_docker = oncotator_docker,
                    gatk_docker = gatk_docker,
                    gatk_override = gatk_override
            }
        }

        # prepare a vcf for evaluation
        File eval_vcf = select_first([MafToVcfConversion.output_vcf, variant])
        File eval_vcf_idx = select_first([MafToVcfConversion.output_vcf_idx, variant_index])

        call Concordance {
            input:
                gatk_override = gatk_override,
                intervals = intervals,
                gatk_docker = gatk_docker,
                truth = MakeTruth.output_vcf,
                truth_idx = MakeTruth.output_vcf_idx,
                eval = eval_vcf,
                eval_idx = eval_vcf_idx,
                prefix = prefix,
                disk_space = truth_vcf_size + disk_padding
        }

        call ConvertToTable as CreateSensitivityTable {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                input_vcf = Concordance.tpfn,
                input_vcf_idx = Concordance.tpfn_idx,
                is_maf_input = is_maf_input,
                disk_space = ceil(size(Concordance.tpfn, "GB") + size(Concordance.tpfn_idx, "GB")) + disk_padding
        }

        call ConvertToTable as CreatePrecisionTable {
            input:
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                input_vcf = Concordance.tpfp,
                input_vcf_idx = Concordance.tpfp_idx,
                bam = bam,
                bam_index = bam_index,
                is_maf_input = is_maf_input,
                disk_space = ceil(size(Concordance.tpfp, "GB") + size(Concordance.tpfp_idx, "GB") + size(bam, "GB") + size(bam_index, "GB")) + disk_padding
        }
    }
    

    call CombineTables as SensitivityTables {
        input:
            input_tables = CreateSensitivityTable.table,
            prefix = prefix + "_sensitivity",
            gatk_docker = gatk_docker,
            disk_space = disk_padding
    }

    call CombineTables as PrecisionTables {
        input:
            input_tables = CreatePrecisionTable.table,
            prefix = prefix + "_precision",
            gatk_docker = gatk_docker,
            disk_space = disk_padding
    }

    call AnalyzeSensitivity {
        input:
            input_table = SensitivityTables.table,
            python_script = python_script,
            prefix = prefix,
            depth_bins = depth_bins,
            depth_bin_width = depth_bin_width,
            python_docker = python_docker,
            disk_space = ceil(size(SensitivityTables.table, "GB")) + disk_padding
    }

    call AnalyzePerformance {
        input:
            precision_table = PrecisionTables.table,
            sensitivity_table = SensitivityTables.table,
            r_script = r_script,
            prefix = prefix,
            depth_bins = depth_bins,
            depth_bin_width = depth_bin_width,
            max_depth = max_depth,
            disk_space = ceil(size(PrecisionTables.table, "GB") + size(SensitivityTables.table, "GB")) + disk_padding
    }

    call CombineTables as SummaryTables {
        input:
            input_tables = Concordance.summary,
            prefix = prefix,
            gatk_docker = gatk_docker,
            disk_space = disk_padding
    }

    call Jaccard as JaccardSNP {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            calls = eval_vcf,
            calls_idx = eval_vcf_idx,
            prefix = prefix,
            type = "SNP"
    }

    call Jaccard as JaccardINDEL {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            calls = eval_vcf,
            calls_idx = eval_vcf_idx,
            prefix = prefix,
            type = "INDEL"
    }

    output {
        File truth_stats = RestrictIntervals.output_stats

        File snp_table = AnalyzeSensitivity.snp_table
        File snp_plot = AnalyzeSensitivity.snp_plot
        File indel_table = AnalyzeSensitivity.indel_table
        File indel_plot = AnalyzeSensitivity.indel_plot

        File expected_estimated_af_plot = AnalyzePerformance.expected_estimated_af_plot
        File fp_estimated_af_plot = AnalyzePerformance.fp_estimated_af_plot
        File depth_roc_plot = AnalyzePerformance.depth_roc_plot

        File summary = SummaryTables.table
        File raw_sensitivity_table = SensitivityTables.table
        File raw_precision_table = PrecisionTables.table

        File snp_jaccard_table = JaccardSNP.table
        File indel_jaccard_table = JaccardINDEL.table
        Array[File] filter_analysis = Concordance.filter_analysis
        Array[File] tpfn = Concordance.tpfn
        Array[File] tpfn_idx = Concordance.tpfn_idx
        Array[File] tpfp = Concordance.tpfp
        Array[File] tpfp_idx = Concordance.tpfp_idx
        Array[File] ftnfn = Concordance.ftnfn
        Array[File] ftnfn_idx = Concordance.ftnfn_idx
    }
} #end of workflow

#### Tasks for making truth
task RestrictIntervals {
    File? gatk_override
    String gatk_docker
    File vcf
    File vcf_idx
    File? intervals

    File genotype_filter_python
    Float? genotype_rate
    Int? genotype_quality_cutoff
    Int? read_depth_cutoff

    Float GT_rate = select_first([genotype_rate, 0])
    Int GQ_cutoff = select_first([genotype_quality_cutoff, 0])
    Int DP_cutoff = select_first([read_depth_cutoff, 0])

    String basename = basename(vcf, ".vcf")

    Int? mem
    Int? disk_space

    Int mem_size = select_first([mem, 4])

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${mem_size}g" SelectVariants \
            -V ${vcf} \
            -O restricted.vcf \
            ${"-L " + intervals}

        python ${genotype_filter_python} \
            --GT-rate ${GT_rate} \
            --GQ-cutoff ${GQ_cutoff} \
            --DP-cutoff ${DP_cutoff} \
            restricted.vcf \
            "${basename}_restricted.vcf"

        gatk --java-options "-Xmx${mem_size}g" IndexFeatureFile -I "${basename}_restricted.vcf"
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: 2
    }

    output {
        File output_vcf = "${basename}_restricted.vcf"
        File output_vcf_idx = "${basename}_restricted.vcf.idx"
        File output_stats = "${basename}_restricted.vcf.stats.txt"
    }
}


task MakeTruth {
    File? gatk_override
    String gatk_docker
    File vcf
    File vcf_idx
    File bam
    File bam_idx
    Int  max_depth

    Int? mem
    Int? disk_space
    
    Int mem_size = select_first([mem, 4])

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        
        gatk --java-options "-Xmx${mem_size}g" CalculateMixingFractions -V ${vcf} -I ${bam} -O mixing.table
        gatk --java-options "-Xmx${mem_size}g" AnnotateVcfWithExpectedAlleleFraction -V ${vcf} -O af_exp.vcf --mixing-fractions  mixing.table
        gatk --java-options "-Xmx${mem_size}g" AnnotateVcfWithBamDepth -V af_exp.vcf -I ${bam} -O bam_depth.vcf
        gatk --java-options "-Xmx${mem_size}g" SelectVariants -V bam_depth.vcf --select "BAM_DEPTH < ${max_depth}" -O truth.vcf
    }

    runtime {
        docker: "${gatk_docker}"
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: 2
    }
    
    output { 
    	File mixing = "mixing.table"
    	File output_vcf = "truth.vcf"
        File output_vcf_idx = "truth.vcf.idx"
    }
}

### Tasks for analysing sensitivity

task ConvertToTable {
  File? gatk_override
  String gatk_docker
  File input_vcf
  File input_vcf_idx
  Boolean is_maf_input
  File? bam
  File? bam_index

  Int? mem
  Int? disk_space

  Int mem_size = select_first([mem, 4])

  command <<<
      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      # Compute TP & FP table
      if [[ ${input_vcf} == *-tpfp.vcf ]]; then

          # Attach BAM_DEPTH information
          gatk --java-options "-Xmx${mem_size}g" AnnotateVcfWithBamDepth -V ${input_vcf} -I ${bam} -O bam_depth_tpfp.vcf

          if [[ ${is_maf_input} == "true" ]]; then

              gatk --java-options "-Xmx${mem_size}g" VariantsToTable \
                   -V bam_depth_tpfp.vcf -F CHROM -F POS -F REF -F ALT -F STATUS -F BAM_DEPTH -F TYPE -GF ref_count -GF alt_count -O tmp.table
              tail -n +2 tmp.table | \
              awk 'BEGIN{OFS="\t"; FS="\t"; print "LOCUS\tSTATUS\tBAM_DEPTH\tAF_EST\tTYPE"}{print $1"_"$2"_"$3"_"$4,$5,$6,$9/($8+$9),$7}' > result.table
          else
              gatk --java-options "-Xmx${mem_size}g" VariantsToTable \
                   -V bam_depth_tpfp.vcf -F CHROM -F POS -F REF -F ALT -F STATUS -F BAM_DEPTH -F TYPE -GF AD -O tmp.table
              tail -n +2 tmp.table | \
              awk 'BEGIN{OFS="\t"; FS="\t"; print "LOCUS\tSTATUS\tBAM_DEPTH\tAF_EST\tTYPE"}{split($NF,c,","); d=c[1]+c[2]; print $1"_"$2"_"$3"_"$4,$5,$6,c[2]/d,$7}' OFS="\t" > result.table
          fi
      # Compute TP & FN table
      else
          gatk --java-options "-Xmx${mem_size}g" VariantsToTable \
               -V ${input_vcf} -F CHROM -F POS -F REF -F ALT -F STATUS -F BAM_DEPTH -F AF_EXP -F TYPE -O tmp.table

          tail -n +2 tmp.table | \
          awk 'BEGIN{OFS="\t"; FS="\t"; print "LOCUS\tSTATUS\tBAM_DEPTH\tAF_EXP\tTYPE";}{print $1"_"$2"_"$3"_"$4,$5,$6,$7,$8}' > result.table
      fi
  >>>

  runtime {
      docker: gatk_docker
      bootDiskSizeGb: 12
      memory: mem_size + " GB"
      disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
      preemptible: 2
  }

  output {
      File table = "result.table"
  }
}

task CombineTables {
    Array[File] input_tables
    String prefix
    String gatk_docker
    Int? disk_space

    command {
        for file in ${sep=' ' input_tables}; do
           head -n 1 $file > header
           tail -n +2 $file >> body
        done

        cat header body > "${prefix}_combined.txt"
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: "1 GB"
        disks: "local-disk " + select_first([disk_space, 20]) + " HDD"
        preemptible: 2
    }

    output {
        File table = "${prefix}_combined.txt"
    }
}

task AnalyzeSensitivity {
    File input_table
    File python_script
    String prefix
    Array[Int] depth_bins
    Int depth_bin_width
    String python_docker
    Int? disk_space
    Int? memory

    command {
        python ${python_script} \
               --prefix ${prefix} --input_file ${input_table} \
               --depth_bins ${sep = ' ' depth_bins} \
               --depth_bin_width ${depth_bin_width}
    }

    runtime {
        docker: python_docker
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        memory: select_first([memory, 8]) + " GB"
        continueOnReturnCode: [0,1]
        preemptible: 2
    }

    output {
        File snp_table = "${prefix}_SNP_sensitivity.tsv"
        File snp_plot = "${prefix}_SNP_sensitivity.pdf"
        File indel_table = "${prefix}_Indel_sensitivity.tsv"
        File indel_plot = "${prefix}_Indel_sensitivity.pdf"
    }
}

task AnalyzePerformance {
    File precision_table
    File sensitivity_table
    File r_script
    String prefix
    Array[Int] depth_bins
    Int depth_bin_width
    Int max_depth

    Int? disk_space
    Int? memory

    command {
        Rscript --vanilla ${r_script} ${prefix} ${sensitivity_table} ${precision_table} ${sep=',' depth_bins} ${depth_bin_width} ${max_depth}
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/tag-tools:0.2.4"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        memory: select_first([memory, 8]) + " GB"
        preemptible: 3
    }

    output {
        File expected_estimated_af_plot = "${prefix}_expected_estimated_af.pdf"
        File fp_estimated_af_plot = "${prefix}_fp_estimated_af.pdf"
        File depth_roc_plot = "${prefix}_depth_roc.pdf"
    }
}



# Make Jaccard index table for SNVs or indels from an array of called vcfs
task Jaccard {
    File? gatk_override
    String gatk_docker
    Array[File] calls
    Array[File] calls_idx
    String prefix
    String type #SNP or INDEL

    Int? mem
    Int? disk_space

    Int mem_size = select_first([mem, 4])

    command <<<
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        result="${prefix}_${type}_jaccard.txt"
        touch $result

        # instead of making soft links, move vcf and the index files
        mkdir vcfs
        mv ${sep = ' ' calls} vcfs/
        VCF_PATH=`ls vcfs/*`
        mv ${sep = ' ' calls_idx} vcfs/

        count=0
        for vcf in $VCF_PATH; do
            ((count++))
            gatk --java-options "-Xmx${mem_size}g" SelectVariants -V $vcf --select-type-to-include ${type} -O ${type}_only_$count.vcf
        done

        for file1 in ${type}_only*.vcf; do
            column=0
            for file2 in ${type}_only*.vcf; do
            ((column++))
            if [ $column != 1 ]; then
                printf "\t" >> $result
            fi

            if [ $file1 == $file2 ]; then
                printf 1.0000 >> $result
            else
                gatk --java-options "-Xmx${mem_size}g" SelectVariants -V $file1 --concordance $file2 -O overlap.vcf
                overlap=`grep -v '#' overlap.vcf | wc -l`

                num1=`grep -v '#' $file1 | wc -l`
                num2=`grep -v '#' $file2 | wc -l`
                just1=$((num1 - overlap))
                just2=$((num2 - overlap))

                total=$((overlap + just1 + just2))
                jaccard=`echo $overlap $total | awk '{print $1 / $2}'`

                printf "%0.6f" $jaccard >> $result
            fi
            done
            printf "\n" >> $result
        done
    >>>

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: mem_size + " GB"
        disks: "local-disk " + select_first([disk_space, 50]) + " HDD"
        preemptible: 2
    }

    output { 
    	File table = "${prefix}_${type}_jaccard.txt" 
    }
}

task Concordance {
      File? gatk_override
      String gatk_docker
      File? intervals
      File truth
      File truth_idx
      File eval
      File eval_idx
      String prefix

      Int? mem
      Int? disk_space

      Int mem_size = select_first([mem, 4])

      Boolean vcf_not_compressed = basename(eval, ".vcf.gz") == basename(eval)
      String vcf_extension = if vcf_not_compressed then "vcf" else "vcf.gz"
      String vcf_index_extension = if vcf_not_compressed then "idx" else "tbi"

      command {
          set -e
          export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

          ln -s ${eval} "input.${vcf_extension}"
          ln -s ${eval_idx} "input.${vcf_extension}.${vcf_index_extension}"

          gatk --java-options "-Xmx${mem_size}g" Concordance ${"-L " + intervals} \
            -truth ${truth} -eval "input.${vcf_extension}" \
            -tpfn "${prefix}-tpfn.vcf" \
            -tpfp "${prefix}-tpfp.vcf" \
            -ftnfn "${prefix}-ftnfn.vcf" \
            --filter-analysis "${prefix}-filter-analysis.txt" \
            -summary "${prefix}-summary.tsv"
      }

      runtime {
          docker: gatk_docker
          bootDiskSizeGb: 12
          memory: mem_size + " GB"
          disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
          preemptible: 2
      }

      output {
            File tpfn = "${prefix}-tpfn.vcf"
            File tpfn_idx = "${prefix}-tpfn.vcf.idx"
            File tpfp = "${prefix}-tpfp.vcf"
            File tpfp_idx = "${prefix}-tpfp.vcf.idx"
            File ftnfn = "${prefix}-ftnfn.vcf"
            File ftnfn_idx = "${prefix}-ftnfn.vcf.idx"
            File filter_analysis = "${prefix}-filter-analysis.txt"
            File summary = "${prefix}-summary.tsv"
      }
}
