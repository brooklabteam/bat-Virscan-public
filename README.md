# bat-Virscan-public

Here is the [README file](https://github.com/brooklabteam/bat-Virscan-public/blob/main/prep-scripts/counts-to-z-score.md) for the pipeline for binning peptide counts off the PhIP-Seq pipeline into Z-scores.
- the actual pipeline needs to be added in this subfolder


Below is the previous README from Emily's folder -- need to edit here.

# bat-VirScan

Here are data and previous scripts for PhIP-Seq runs of Singaporean *Eonycteris spelaea* bats (raised in captivity; only 5 bats) and Australian *Pteropus alecto* bats (wild-caught but rehabilitated; 77 bats).

The files available are organized as:

- Two raw count data prior to AVARDA can be found in "AVARDA" > "input"
- Two folders with AVARDA output for each dataset, named "Espelaea_zscores_binning_avarda 09-December-2021_11-00" and "Palecto_zscores_binning_avarda 09-December-2021_11-11".
- Two z-score files (post AVARDA) for each of the PhIP-Seq runs (one per species): "Espelaea_clean_results_summary.csv" and "Palecto_clean_results_summary.csv"
- The "All_bat_sum9Dec2021.csv" file that has the combined P.alecto and E.spelaea VirScan data with metadata
- One metadata file (both species): "bat_metadata_9Dec2021.csv"


Here is a brief overview of how we obtained viral hit data:
1. Raw sequences were aligned using the Tuxedo suite and the information found in the "Raw Sequence alignment" > "Scripts" folder.
2. Raw count data were then input into the "counts to Z-scores_Brooklab.Rmd" file in order to obtain z-scores from a ranked binning approach (Mina et al. 2019).
3. Final adjusted z-scores were then input into the AVARDA pipeline. Please see the Readme file for a description.
4. The AVARDA output probablistic viral hits and those results can be found in "All_bat_sum9Dec2021.csv".
