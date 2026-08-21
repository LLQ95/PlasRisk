@echo off
echo.
echo PlasRisk has been installed successfully.
echo.
echo For full annotation capability (ARG, VF, replicon, biocide/metal),
echo install abricate and blast:
echo.
echo     conda install -c bioconda -c conda-forge abricate blast
echo.
echo Then set up abricate databases:
echo.
echo     abricate-get_db --db card --force
echo     abricate-get_db --db vfdb --force
echo     abricate-get_db --db plasmidfinder --force
echo.
echo Quick start:
echo     plasrisk plasmid.fasta
echo.
