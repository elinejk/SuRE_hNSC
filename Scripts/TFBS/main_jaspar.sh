cd /gpfs/home6/evzanten/jaspar_enrichment/scripts_enrichment/

for tfbs in YY1 GABPA TP53 TP63 ZBTB33 ; do sbatch "/gpfs/home6/evzanten/jaspar_enrichment/scripts_enrichment/LOLA_JASPAR_twoSets.sh" $tfbs ; done