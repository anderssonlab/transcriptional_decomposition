# transcriptional_decomposition

This README is a guide to the R scripts included with our manuscript entitled "Transcriptional decomposition reveals active chromatin architectures and cell specific regulatory interactions". Further questions regarding code should be made to sarah@binf.ku.dk or robin@binf.ku.dk.

1. transcriptional_decomposition.R - script for getting from raw CTSS CAGE samples to decomposed PD/PI components

To generate transcriptionally decomposed components (PD/PI or often referred to as RW/IID), the following information is required:

* Raw RNA tag counts across the genome (or chromosomes of interest) in the form of CTSS files (CAGE transcription start site), or (non-sparse) count data which can be aggregated into fixed width bin quantities.

* A chosen bin width (10kb recommended for human, larger for sparse data and less for compact genomes).

* Library depths of samples.

The script first compiles the individual samples into a genomic range (using GRanges) containing all bps with non-zero counts, which is used to calculate an aggregated count across fixed bin widths over the given chromosomes/genome. 

To run the transcriptional decomposition, R-INLA (Rue et al 2009) is required. The model assumes a negative binomial distribution on the dataset and that the log of the counts are composed of two random effects, a random walk and an independent component. The R-INLA call is made in two steps; since the prior is vague and the data often sparse (i.e. at centromeres etc.), the model may get stuck, so a first run (adding 1 to the diagonal of the variance-covariance matrix) helps to find more stable priors, and allows the second run to more efficiently generate the correct results.

The R-INLA call saves out the results to the given path; these are reopened and the RW and IID components extracted for each cell type.

2. generate_xad_boundaries.R - script for calling XAD boundaries, calls out to getRankedBounds (3.)

3. getRankedBounds.R - function for generating XAD boundaries from RW decomposed vectors

Implements the boundary finding algorithm (see methods in the paper, and/or illustration in fig. 4A in the main paper)

4. generate_interactions_dataset.R - function for generating annotated bin-bin pairs for interaction modelling

Generates a dataframe of all intrachrosomal bin-bin pairs within (default) 2MB, each with annotated expression parameters (see paper supplementary for list) and interaction significance scores measured for GM12878 CaptureHIC (Mifsud et al 2015) using CHiCAGO (Cairns et al 2016).

5. predict_interactions_on_ENCODE.R - script which takes annotated bin-bin pairs, models against GM12878 interaction counts via a random forest, and predicts on other ENCODE cell lines

6. [not added yet] DE_Hela_GM12878.R - transcriptional decomposition script to include posterior estimates of the difference between two cell lines
