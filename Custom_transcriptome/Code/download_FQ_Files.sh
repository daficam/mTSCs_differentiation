

https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=PRJEB9375&result=read_run&fields=study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,sra_ftp,sra_galaxy,cram_index_ftp,cram_index_galaxy&download=txt



cut -f 12 PRJEB9375.txt | grep -v submitted | sed 's/^/wget --ftp-user=anonymous --continue --no-host-directories /g' > wgets.sh

chmod 755 wgets.sh

./wgets.sh &> wgets.log &
