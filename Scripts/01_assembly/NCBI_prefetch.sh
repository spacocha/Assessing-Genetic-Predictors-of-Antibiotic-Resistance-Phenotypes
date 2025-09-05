#Obtain isolates from NCBI under BioProject PRJNA664939 (search for BioProject, click on associated 308 SRA entries -> Send Results to Run Selector
#Search for *Pseudomonas* = 220 different SRAs-> Return metadata for this list
#Compare manually with the results from Table S6 from Hannah’s paper.
#All of the ones on the list are present, some extra SRAs that aren’t on the list are available.
#Ones without phenotypes are removed
#This results in 207 isolates from the list of phenotype that is also Pseudomonas.
ml sra-tools/3.0.3
mkdir PREFETCH_dir
prefetch --option-file Accession_list_207_NCBI_isolates.txt --output-directory PREFETCH_dir > prefetch_log.txt
fasterq-dump --split-files PREFETCH_dir/SRR12761*/*.sra

#This results in all of the fastq files for each isolate

