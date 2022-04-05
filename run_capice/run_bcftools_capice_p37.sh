#!/bin/bash

capice_bcftools() {
  local -r header="%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%SYMBOL\t%SYMBOL_SOURCE\t%Gene\t%Feature\t%Feature_type\t%cDNA_position\t%CDS_position\t%Protein_position\t%Amino_acids\t%STRAND\t%SIFT\t%PolyPhen\t%EXON\t%INTRON\t%SpliceAI_pred_DP_AG\t%SpliceAI_pred_DP_AL\t%SpliceAI_pred_DP_DG\t%SpliceAI_pred_DP_DL\t%SpliceAI_pred_DS_AG\t%SpliceAI_pred_DS_AL\t%SpliceAI_pred_DS_DG\t%SpliceAI_pred_DS_DL"
  local -r capiceInputPathHeaderless="capice_input.tsv.headerless"

  local args=()
  args+=("+split-vep")
  args+=("-d")
  args+=("-f" "${header}\n")
  args+=("-o" "${capiceInputPathHeaderless}")
  args+=("lp37.vcf.gz")

  singularity exec /groups/solve-rd/tmp10/mslofstra/vip/images/bcftools-1.14.sif bcftools "${args[@]}"

  echo -e "${header}" | cat - "${capiceInputPathHeaderless}" > "capice_input_p37.tsv"
}

export SINGULARITY_BIND=/apps,/groups,/tmp
capice_bcftools