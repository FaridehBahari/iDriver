#!/bin/bash
 
# List of cohorts

cohorts=( #"Pan_Cancer"
 "Pancan-no-skin-melanoma-lymph"
"Biliary-AdenoCA"  "Bladder-TCC" 
"Bone-Leiomyo"  "Bone-Osteosarc"
"Breast-AdenoCa"    "CNS-GBM"
"CNS-Medullo" "CNS-Oligo"  "CNS-PiloAstro"    "ColoRect-AdenoCA"
"Eso-AdenoCa"  "Head-SCC"   "Kidney-RCC" "Liver-HCC"
 "Lung-AdenoCA" "Lung-SCC" "Lymph-BNHL"
"Lymph-CLL"        "Myeloid-MPN"      "Ovary-AdenoCA"
"Panc-AdenoCA"     "Panc-Endocrine"   "Prost-AdenoCA"    "Skin-Melanoma"
"Stomach-AdenoCA"       "Uterus-AdenoCA" "Kidney-ChRCC"
 "Cervix-SCC" "Thy-AdenoCA"
"Breast-DCIS"  "Breast-LobularCa"   "Lymph-NOS"
)

# cohorts=( 
  # "Pancan-no-skin-melanoma-lymph" "ColoRect-AdenoCA" "CNS-GBM"
  # "Head-SCC"         "Skin-Melanoma"         "Lung-SCC" "Lung-AdenoCA"
  # "Stomach-AdenoCA"  "Uterus-AdenoCA"   "Biliary-AdenoCA"  "Eso-AdenoCa" "Lymph-BNHL"
# )

###### Rho = 0 ###### 
# Create "../extdata/procInput/iDriverInputs/Pan_CancerRho0FILE.tsv" by replacing all Rhos from Pan_Cancer cohort of withHypers to 0
template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 40 "{cohort}" "../extdata/output/observed/Rho0/rmHyperMuts/pooled_betaRho0_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/procInput/iDriverInputs/Pan_CancerRho0FILE.tsv"'
# template_command='Rscript main/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/Rho0/betaRho0_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/procInput/iDriverInputs/Pan_CancerRho0FILE.tsv"'

###### population-level ###### 
# Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/population_level/withHyperMuts/pooled_betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "population_level" 
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/population_level/rmHyperMuts/pooled_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "population_level" "../extdata/output/observed/population_level/withHyperMuts/pooled_betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript main/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/simulated/population_level_rmHyperMuts/betaRhoPancan_ElemType/" {flag} TRUE 4 "population_level" "../extdata/output_release2.0/observed/population_level/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### Burden-free ###### 
# Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Burden-free"
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/Burden_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Burden-free" "../extdata/output/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript main/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/simulated/Burden_free/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Burden-free" "../extdata/output_release2.0/observed/Burden_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### FI-free (version 1) ######
# Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "FI-free"
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 50 "{cohort}" "../extdata/output/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "FI-free" "../extdata/output/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript main/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/FI_free/rmHyperMuts_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "FI-free" "../extdata/output_release2.0/observed/FI_free/withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'



###### Original (withHyperMuts) ######
# Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/Pan_Cancer/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original"
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" {flag} FALSE 4 "Original" "../extdata/output/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'



###### Original (rmHyperMuts) ######
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'
# template_command='Rscript main/iDriver.R "simulated" "../extdata/procInput/iDriverInputs/simulated/Sanger/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/simulated/Sanger_rmHyperMuts/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/simulated/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/simulated/rmHyperMuts_betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###### Original (rmHyperMuts_GBM pElems) ######
# Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed_withHyperMutated/Pan_Cancer/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output_release2.0/observed/pooled_withHyperMuts/gbm_pElems/" FALSE FALSE 4 "Original"
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 75 "{cohort}" "../extdata/output_release2.0/observed/pooled/gbm_pElems/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/gbm_pElems/observed_Pan_Cancer.tsv"'


# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs_2024/observed/{cohort}/pElems/eMET_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output_release2.0/observed/pooled/forAverageTis/" {flag} TRUE 4 "Original" "../extdata/output_release2.0/observed/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv"'

###################### Ablation (BMR effect) ######################################3
# DP
# nohup Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/DP/withHyperMutated/Pan_Cancer/pElems/DP_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/ablation/BMReffect/DP/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original" > ../ablation_DPbmr &
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/DP/rmHyperMutated/{cohort}/pElems/DP_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/ablation/BMReffect/DP/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output/observed/ablation/BMReffect/DP/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv" '

# Baseline
# nohup Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/baseLine_localRates/withHyperMutated/pElems/Pan_Cancer_Baseline_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/ablation/BMReffect/baseLine_localRates/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original" > ../ablation_BaselineBMR &
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/baseLine_localRates/rmHyperMutated/pElems/{cohort}_Baseline_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/ablation/BMReffect/baseLine_localRates/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output/observed/ablation/BMReffect/baseLine_localRates/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv" '

# GBM
# nohup Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/withHyperMutated/Pan_Cancer/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/ablation/BMReffect/GBM/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original" > ../ablation_GBMbmr &
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/pElems/GBM_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/ablation/BMReffect/GBM/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output/observed/ablation/BMReffect/GBM/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv" '

# Dig
# nohup Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/Dig/withHyperMutated/Pan_Cancer/pElems/Dig_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/Pan_Cancer_CADD6_elemTypeParams.tsv" "Pan_Cancer" 70 "Pan_Cancer" "../extdata/output/observed/ablation/BMReffect/Dig/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/" FALSE FALSE 4 "Original" > ../ablation_DigBMR &
# template_command='Rscript main/iDriver.R "observed" "../extdata/procInput/iDriverInputs/cadd_scores.tsv" "../extdata/procInput/BMRs/observed/Dig/rmHyperMutated/{cohort}/pElems/Dig_orig_pElmens.tsv" "../extdata/procInput/Mu_sd/elemType/observed/{cohort}_CADD6_elemTypeParams.tsv" "{cohort}" 70 "{cohort}" "../extdata/output/observed/ablation/BMReffect/Dig/pooled/betaRhoPancan_lengthCat_quartilesElemType/" {flag} TRUE 4 "Original" "../extdata/output/observed/ablation/BMReffect/Dig/pooled_withHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pan_Cancer.tsv" '





# Loop over cohorts
for cohort in "${cohorts[@]}"
do
  echo "Running analysis for cohort: $cohort"
 
  # Determine logical flag for excluding Lymphoma and Melanoma
  if [[ "$cohort" == "Pancan-no-skin-melanoma-lymph" ]]; then
    second_last_flag="TRUE"
  else
    second_last_flag="FALSE"
  fi
 
  # Replace placeholders in the template
  modified_command="${template_command//\{cohort\}/$cohort}"
  modified_command="${modified_command//\{flag\}/$second_last_flag}"
 
  # Execute the command
  eval $modified_command
done