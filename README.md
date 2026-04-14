# 2025-prediction-contest

The objective will be to use data available on The Triticeae Toolbox [T3/Wheat](https://wheat.triticeaetoolbox.org/) to predict yield (Grain yield - kg/ha|CO_321:0001218) performance of wheat accessions in 9 separate test trials across the USA. Important metadata for these trials is on T3/Wheat, but the phenotypes measured in them are not. For all of these trials, genotypic data is available, though it is somewhat messy: 1. Not all accessions in the trial are genotyped; 2. Different genotyping protocols have been used in different trials.

The names of the 9 trials are:

    AWY1_DVPWA_2024
    TCAP_2025_MANKS
    25_Big6_SVREC_SVREC
    OHRWW_2025_SPO
    CornellMaster_2025_McGowan
    24Crk_AY2-3
    2025_AYT_Aurora
    YT_Urb_25
    STP1_2025_MCG

These trials are also available as the public list named T3 Prediction Challenge Trials on T3/Wheat.

Contestants will submit a folder containing predictions for genotyped entries in the trials for "CV0" and "CV00" cases (Jarquin et al. 2017; see below). The winning algorithm will be judged on the basis of the highest average prediction accuracy across trials and cross validation cases.
Training Data

We have identified T3 data that could be reasonably applied to predicting the Predictathon trials using simple heuristics:

    We identified phenotyping trials that evaluated some of the same accessions as each Predictathon trial. The assumption is that the overall germplasm in those trials will be related to the Predictathon accessions.
    We identified archived VCFs with genotype data from a minimum number of each of the accessions in the putative training trials.

This data was generated from scripts found in this GitHub repository:
https://github.com/jeanlucj/T3_predictathon_find_training_trials
~                                                                      
