#!/bin/bash

capice_vep() {
  local args=()
  args+=("--input_file" "/groups/solve-rd/tmp10/mslofstra/vip/test/resources/lp.vcf.gz")
  args+=("--format" "vcf")
  args+=("--output_file" "/groups/solve-rd/tmp10/mslofstra/lp37.vcf.gz")
  args+=("--vcf")
  args+=("--compress_output" "bgzip")
  args+=("--no_stats")
  args+=("--fasta" "/groups/solve-rd/tmp10/projects/vip/git/vip/resources/GRCh37/human_g1k_v37.fasta.gz")
  args+=("--offline")
  args+=("--cache")
  args+=("--dir_cache" "/groups/solve-rd/tmp10/mslofstra/vip/resources/vep/cache")
  args+=("--species" "homo_sapiens")
  args+=("--assembly" "GRCh37")
  args+=("--refseq")
  args+=("--exclude_predicted")
  args+=("--use_given_ref")
  args+=("--symbol")
  args+=("--flag_pick_allele")
  args+=("--sift" "s")
  args+=("--polyphen" "s")
  args+=("--total_length")
  args+=("--shift_3prime" "1")
  args+=("--allele_number")
  args+=("--numbers")
  args+=("--dont_skip")
  args+=("--allow_non_variant")
  args+=("--buffer_size" "1000")
  args+=("--fork" "8")
  args+=("--dir_plugins" "/groups/solve-rd/tmp10/mslofstra/vip/resources/vep/plugins")
  args+=("--plugin" "SpliceAI,snv=/groups/solve-rd/tmp10/projects/vip/git/vip/resources/GRCh37/spliceai_scores.masked.snv.hg19.vcf.gz,indel=/groups/solve-rd/tmp10/projects/vip/git/vip/resources/GRCh37/spliceai_scores.masked.indel.hg19.vcf.gz")

  singularity exec /groups/solve-rd/tmp10/mslofstra/vip/images/vep-105.0.sif vep "${args[@]}"
}

export SINGULARITY_BIND=/apps,/groups,/tmp
capice_vep