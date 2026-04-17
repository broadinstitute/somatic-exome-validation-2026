# This workflow only works for converting hg19 MAF to hg19 VCF.
# This code is provided by Lee Lichtenstein <lichtens@broadinstitute.org>

workflow MafToVcfConversion {
    File input_maf
    String sample_id
    
    # python script to extract major columns in MAF to avoid errors
    File petite_maf_python

    # oncotator aux files
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? oncotator_exe
    File? default_config_file

    # docker files
    String? oncotator_docker
    String? gatk_docker
    File? gatk_override

    call MAFtoVCF {
        input:
            input_maf = input_maf,
            sample_id = sample_id,
            petite_maf_python = petite_maf_python,
            onco_ds_tar_gz = onco_ds_tar_gz,
            onco_ds_local_db_dir = onco_ds_local_db_dir,
            oncotator_exe = oncotator_exe,
            default_config_file = default_config_file
    }

    call GenerateVcfIndex {
        input:
           input_vcf = MAFtoVCF.output_vcf,
           gatk_docker = gatk_docker,
           gatk_override = gatk_override
    }

    output {
        File output_vcf = MAFtoVCF.output_vcf
        File output_vcf_idx = GenerateVcfIndex.output_vcf_idx
    }
}

task MAFtoVCF {
    File input_maf
    String sample_id
    File petite_maf_python

    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? oncotator_exe
    File? default_config_file
    
    # runtime
    String? oncotator_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    command <<<
        set -e

        # local db dir is a directory and has been specified
        if [[ -d "${onco_ds_local_db_dir}" ]]; then
            echo "Using local db-dir: ${onco_ds_local_db_dir}"
            echo "THIS ONLY WORKS WITHOUT DOCKER!"
            ln -s ${onco_ds_local_db_dir} onco_dbdir
        elif [[ "${onco_ds_tar_gz}" == *.tar.gz ]]; then
            echo "Using given tar file: ${onco_ds_tar_gz}"
            mkdir onco_dbdir
            tar zxvf ${onco_ds_tar_gz} -C onco_dbdir --strip-components 1
        else
            echo "Downloading and installing oncotator datasources from Broad FTP site..."
            # Download and untar the db-dir
            wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz
            tar zxvf oncotator_v1_ds_April052016.tar.gz
            ln -s oncotator_v1_ds_April052016 onco_dbdir
        fi
        
        TARGET_MAF="`basename ${input_maf} .maf`.maf"
        python ${petite_maf_python} ${input_maf} $TARGET_MAF

        ${default="/root/oncotator_venv/bin/oncotator" oncotator_exe} --infer_genotypes true --db-dir onco_dbdir/ $TARGET_MAF ${sample_id}.vcf -o VCF hg19
            ${"--default_config " + default_config_file}
    >>>

    runtime {
        docker: select_first([oncotator_docker, "broadinstitute/oncotator:1.9.6.1"])
        memory: select_first([mem, 3]) + " GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
    }

    output {
        File output_vcf = "${sample_id}.vcf"
    }
}

task GenerateVcfIndex {
    File input_vcf
    String? gatk_docker
    File? gatk_override

    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    
    String vcf_basename = basename(input_vcf)
    
    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        mv ${input_vcf} .

        gatk --java-options "-Xmx3g" IndexFeatureFile -I ${vcf_basename} -O ${vcf_basename}.idx
    >>>

    runtime {
        docker: select_first([gatk_docker, "us.gcr.io/broad-gatk/gatk:4.0.9.0"])
        memory: select_first([mem, 3]) + " GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + select_first([disk_space, 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }
    output {
        File output_vcf_idx = "${vcf_basename}.idx"
    }
}
