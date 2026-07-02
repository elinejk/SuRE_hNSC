# General analyses of the data

## 1. Replicate correlations and qqplot
**replicate_correlation.R:** this script is used to determine the correlation for activity values (average activity for all fragments covering an allele) between the two hNSC replicates.

**replicate_correlation_stratified.R:** this script splits SNPs based on their p-value, expression levels and coverage. For each bin, it then determines the correlation for activity values (average activity for all fragments covering an allele), the correlation for effect sizes (as log2(REF/ALT)), and the concordance for effect sizes between the two hNSC replicates.

**qqplot_pvalues.R:** creates a qqplot for measured (hNPC.wilcoxon.pvalue) and expected (hNPC.wilcoxon.pvalue.random) p-values, using 100,000 randomly selected SuRE SNPs (Figure S1A).

## 2. Minor allele frequency
**MAF.sh:** bash script used to determine the MAF for each SuRE SNP

**MAF_analyses.R:** calculates the minimum, maximum, mean, and median MAF for emVars, control SNPs and the full SuRE dataset

## 3. emVar and control SNP characteristics
Various scripts that were used to compare emVars and control SNPs based on their expression, coverage, MAF and p-values

**emvar_vs_controls.R:** compared emVars and control SNPs based on their expression, fragment count, and MAF. It also creates the plots used for Figures S1C-F

**storeys_pi1.R:** determines the storey's pi1 between emVars of the three cell lines.
