# COVID clustering
This folder contains all the data and code used in the manuscript entitled "Transcriptomic clustering of critically ill COVID-19 patients", published in the European Respiratory Journal (https://erj.ersjournals.com/content/early/2022/08/25/13993003.00592-2022). The preprint is available at https://www.medrxiv.org/content/10.1101/2022.03.01.22271576v1.

The folder contains the following files:

- "covid_clustering.R": R code. This file uploads all the following data files, that must be in the same directory.

- "transcript_data.Rdata": Transcriptomic data. FastQ files available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197204.

- "transcript_data_trt.Rdata": Transcriptomic data at ICU day 4. FastQ files available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206264.

- "miRNA_counts.csv": miRNA raw counts. FastQ files available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197258.

- "molecules_downnetwork.txt": Results from Ingenuity Pathway Analysys.

- "molecules_upnetwork.txt": Results from Ingenuity Pathway Analysys.

- "clinical_data.csv": Anonymized clinical data.

- "index_data_trt.csv": File containing cluster and treatment received by patients with a Day 4 transcriptome.

## Contact: 

Follow us: https://twitter.com/crit_lab

Cecilia López-Martínez (clm@crit-lab.org), Guillermo M Albaiceta (gma@crit-lab.org), Laura Amado-Rodríguez (lar@crit-lab.org).
