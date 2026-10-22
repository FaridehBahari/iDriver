import os
import subprocess

# List of cohorts
# 
cohorts = [ 
             "Myeloid-MPN", "Biliary-AdenoCA", 'Pancan-no-skin-melanoma-lymph', "Liver-HCC",
             "Panc-AdenoCA", "Prost-AdenoCA", "Pan_Cancer"
             "Cervix-SCC",  "CNS-Oligo",
             "Bladder-TCC", "Bone-Leiomyo", "Bone-Osteosarc",
             "Breast-AdenoCa", "CNS-GBM", "CNS-Medullo", "CNS-PiloAstro",
             "ColoRect-AdenoCA", "Eso-AdenoCa", "Head-SCC", "Kidney-ChRCC",
             "Kidney-RCC", "Lung-AdenoCA", "Lung-SCC", "Lymph-BNHL",
             "Lymph-CLL", "Ovary-AdenoCA",  "Panc-Endocrine",
             "Skin-Melanoma", "Stomach-AdenoCA", "Thy-AdenoCA",
             "Uterus-AdenoCA",
  
   # 'Pancan-no-skin-melanoma-lymph', "Pan_Cancer", 'Biliary-AdenoCA',
   #  'Skin-Melanoma',  'CNS-GBM', 'ColoRect-AdenoCA',  'Eso-AdenoCa',
   #  'Head-SCC', 'Lung-AdenoCA', 'Lung-SCC',
   #  'Stomach-AdenoCA',  'Uterus-AdenoCA', "Lymph-BNHL"
  
  
]


# Path to the template ini file
template_ini_path = 'configs/sim_setting_finetuning_iDriver2026.ini'

# Read the template ini file
with open(template_ini_path, 'r') as file:
    ini_content = file.read()

# Loop over the cohorts and run the algorithm
for cohort in cohorts:
    print(cohort)
    # Replace the cohort name in the ini content
    modified_ini_content = ini_content.replace('Pancan-no-skin-melanoma-lymph', cohort)
    
    # Save the modified ini content to a new file
    modified_ini_path = f'configs/sim_setting_finetuning_iDriver2026_{cohort}.ini'
    with open(modified_ini_path, 'w') as file:
        file.write(modified_ini_content)
    
    # Run the algorithm with the modified ini file
    command = f'python -u run_eMET.py {modified_ini_path} ../../iDriver/extdata/procInput/BMRs/observed/rmHyperMutated/{cohort}/GBM/GBM_model.pkl 100'
    # command = f'python -u run_eMET.py {modified_ini_path} ../../iDriver/extdata/procInput/BMRs/observed/withHyperMutated/{cohort}/GBM/GBM_model.pkl 100'
    
    subprocess.run(command, shell=True)

# from BMR directory run the following command:
# nohup python -u ../../iDriver/iDriver/prepareInput/6_run_eMET_allCohorts.py > ../../run_allCohorts_eMET &

