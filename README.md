# bat-Virscan-public

## Data preparation

We received data in raw peptide count form from 5 Singaporean *Eonycteris spelaea* bats (raised in captivity) and from 77 Australian *Pteropus alecto* bats (wild-caught but rehabilitated) after Phage Immunoprecipitation Sequencing (PhIP-Seq) was run in the Wang lab - see subfolders within "/raw-data/PhIP-Seq/" to access these counts. We first sought to convert these peptide counts into Z-scores for input into the AVARDA pipeline (*link to repo*) for rank categorization of overlapping peptides into putative viral exposures. See [here](https://github.com/brooklabteam/bat-Virscan-public/blob/main/prep-scripts/counts-to-z-score.md) for documentation of the process and pipeline for binning peptide counts off the PhIP-Seq pipeline into Z-scores.

Following establishment of these Z-scores, we input these files into the AVARDA pipeline (*link to paper or repo here*) to establish putative viral exposures. See "AVARDA" subfolder within the "raw-data" folder for outputs from this pipeline.

After running AVARDA, we merged these viral exposure data with metadata provided by collaborators for all subsequent analyses. The merged AVARDA output data + metadata for both species of bat is available in the "working-data" subfolder as file "merge-data.csv".

We reference these data in all subsequent figure development files.

## Repo

This repo is organized by figure. Background scripts for generating datasets are found in the "prep-scripts" folder, while all code needed to produce analyses and generate plots for each figure are in the corresponding 'fig1', 'fig2', 'fig3', and 'fig4' subfolders. These folders output files into 'final-figures' and 'supp-figures' accordingly..