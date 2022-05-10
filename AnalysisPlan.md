GWAS of PC1 of white matter hyperintensities and insular cortical
thickness
================
Ariana Tang (<ariana.tang@sickkids.ca>)Jean Shin
(<jean.shin@sickkids.ca>)Tomas Paus (<tomas.paus@umontreal.ca>)Zdenka
Pausova (<jean.shin@sickkids.ca>)

------------------------------------------------------------------------

-   **Aim**

    *To identify genetic loci associated with the shared variance
    between white matter hyperintensities (WMH) and thickness of the
    insular cortex*

    Note: This project is to be conducted in 2 stages: (Stage 1) to
    calculate the first principal component (PC1) scores for white
    matter hyperintensities (WMH) and thickness of the insular cortex,
    and (Stage 2) to perform genome-wide association tests for PC1. For
    Step 1, R scripts will be distributed.

-   **Subject exclusion criteria**

    -   Dementia (≥mild severity) at time of MRI-scanning
    -   Stroke at time of MRI-scanning (use exclusion criteria if
        available for study, either based on clinical data or large
        artery strokes/lacunes in region of interest on MRI)
    -   Multiple sclerosis (if available)
    -   Brain surgery (if available)
    -   Morphological abnormalities (e.g., cysts, brain tumors)
    -   Poor technical quality  
    -   Non-White/non-Caucasian Participants

-   **Brain outcomes**

    -   White matter hyperintensities (WMH): T1, T2 or FLAIR, total load
        as a quantitative variable is preferable, but if not available,
        total load as a semiquantitative variable can be used.
    -   Cortical thickness (CT) averaged between the left and the right
        hemispheres in the 34 FreeSurfur-parcellated regions.
    -   First principal of WMH and thickness of the insular cortex.

-   **Covariates**

    1.  Age: linear and quadratic terms
    2.  Sex
    3.  Intracranial volume (ICV) or brain size
    4.  Current smoking (yes/no)
    5.  Hypertension (yes/no)
    6.  Type 2 diabetes (yes/no)
    7.  Cohort specific covariates: e.g., MR-scanner, etc.
    8.  Genetic principal components (cohort-specific numbers): **Step 2
        only**

------------------------------------------------------------------------

-   **Stage 1 (run with the provided scripts)**

    1.  Descriptive statistics will be calculated for the brain
        phenotypes and covariates (#1-#7)
    2.  Brain outcomes will be first inverse-normal transformed (INT)
        and then adjusted for covariates:
    3.  Models M1 and M2 will be fitted within each sex

    **M1**: INT(WMH) \~ Age + Age<sup>2</sup> + ICV or Brain Volume +
    Hypertension + Type 2 diabetes + Current smoking status +
    Cohort-specific covariates

    **M2**: INT(insulaCT) \~ Age + Age<sup>2</sup> + Hypertension + Type
    2 diabetes + Current smoking status + Cohort-specific covariates

    4.  PC1 scores for WMH and insular CT will be calculated using the
        residuals obtained from fitting M1 and M2.
    5.  Use the PC1 scores returned after running Stage 1: the file is
        in ‘files_for_GWAS’ folder

------------------------------------------------------------------------

-   **Stage 2: Genome-wide association study (in the sex-combined
    participants)**

    1.  Run GWAS under model (M3) using the PC1 scores returned from
        running the Stage 1 scripts (it will be available from the
        outputted ‘GWAS_file/covariate_adjusted_phenotypes_for_GWAS.tsv’
        )

    **M3**: PC1 sim SNP + Age + Age<sup>2</sup> + Sex + ICV or Brain
    Volume + Hypertension + Type 2 diabetes + Current smoking status +
    Cohort-specific covariates + Genetic PCs

\[should I INT(PC1)??\]

2.  Format the GWAS result file with the following column names:

| Column header     | Description                                                     |
|-------------------|-----------------------------------------------------------------|
| ID                | The dbSNP or 1000 Genome Project ID number                      |
| CHROM             | Chromosome number                                               |
| POS               | Base pair position of SNP on the chromosome, Build 37 (or hg19) |
| CODED_ALLELE      | Allele for which the effect has been estimated                  |
| NON_CODED_ALLELE  | Non-coded allele                                                |
| CODED_ALLELE_FREQ | Allele frequency of the coded allele                            |
| N                 | Number of individuals used in the test                          |
| BETA              | Effect size estimate for the coded allele                       |
| SE                | Standard error of the effect size estimate                      |
| P                 | P-value of the effect size                                      |
