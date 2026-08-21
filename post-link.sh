#!/bin/bash
# post-link.sh - Displayed after conda install on Linux/macOS
cat <<EOF

PlasRisk has been installed successfully.

For full annotation capability (ARG, VF, replicon, biocide/metal),
install abricate and blast:

    conda install -c bioconda -c conda-forge abricate blast

Then set up abricate databases:

    abricate-get_db --db card --force
    abricate-get_db --db vfdb --force
    abricate-get_db --db plasmidfinder --force
    abricate-get_db --db resfinder --force
    abricate-get_db --db ncbi --force

Quick start:
    plasrisk plasmid.fasta
    plasrisk *.fasta

For BacMet database setup, see:
    https://github.com/tseemann/abricate#making-your-own-database

EOF
