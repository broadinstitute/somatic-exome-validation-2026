#!/bin/bash
# author: Junko Tsuji <jtsuji@broadinstitute.org>
# Processing read coverage and linking the coverage with GC%

# pan-cancer gene list
GENE_LIST="data/396genes.list"
# overlapping genomic regions between v2 and v6
INTERVALS="data/TWIST-ICE_overlap.bed"
# gencode gtf file version 19 downloaded from https://www.gencodegenes.org/human/release_19.html
GENCODE="gencode.v19.annotation.gtf"

# 1: extract gene intervals
grep -f ${GENE_LIST} -F ${GENCODE} | \
  grep 'tag "basic"' | \
  grep "exon_number 1;" | \
  awk '{if($3 == "exon"){print $0}}' | \
  awk -F "\t" '{split($9,tmp,"gene_name ");
                split(tmp[2],gene,";");
                split($9,tmp,"gene_type ");
                split(tmp[2],gene_type,";");
                print $1,$4-1,$5,gene[1],gene_type[1],$7}' OFS="\t" | \
  grep -v processed_transcript | \
  awk '{print $1,$2,$3,$4":"$6}' OFS="\t" | \
  sed 's/"//g' | sed 's/^chr//g' | \
  bedtools intersect -a stdin -b ${INTERVALS} | \
  sort -k4,4 -k1,1 -k2,2n | \
  python scripts/merge_entries.py - > exon1_pancan_396gene.bed

# 2: convert bed format to interval_list
gatk-4.1.0.0/gatk BedToIntervalList -I exon1_pancan_396gene.bed \
                                    -O exon1_pancan_396gene.interval_list \
                                    -SD /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.dict

# 3: annotate GC% per genomic region
gatk-4.1.0.0/gatk AnnotateIntervals --interval-merging-rule OVERLAPPING_ONLY \
                                    -L exon1_pancan_396gene.interval_list \
                                    -R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
                                    -O exon1_pancan_396gene.gc.txt

# 4: format GC% output into bed file
grep -v ^@ exon1_pancan_396gene.gc.txt | \
  grep -v ^CONTIG | \
  awk '{print $1,$2-1,$3,$4}' OFS="\t" > exon1_pancan_396gene.gc.bed && rm exon1_pancan_396gene.gc.txt

# 5: attach GC% to the first exons of the pancancer genes
bedtools intersect -a exon1_pancan_396gene.bed -b exon1_pancan_396gene.gc.bed -wo | \
  cut -f1-4,8 > tmp && mv tmp exon1_pancan_396gene.bed

# 6: run coverage analysis and attach GC%
for i in `cat *.bam`;
do
  gatk DepthOfCoverage -I $i \
                       -L exon1_pancan_396gene.interval_list \
                       -O `basename $i .bam`.cov \
                       -R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta
  grep -v ^Locus `basename $i .bam`.per_base_coverage.txt | \
  awk '{split($1,pos,":"); $1=pos[1]"\t"pos[2]-1"\t"pos[2]; print $0}' OFS="\t" | \
  bedtools intersect -a stdin -b exon1_pancan_396gene.bed -wo > out
  python scripts/tally_covered_bases.py out > `basename $i .bam`.covered_bases.txt && rm out
done
