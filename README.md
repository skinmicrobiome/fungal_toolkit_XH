# fungal_toolkit_XH
Tools and databases for fungal ITS1 classification and Candida auris clade assignment

Reference: **Huang et al., "Skin metagenomic sequence analysis of Candida auris outbreaks from American healthcare facilities", manuscript in preparation**

**Under Construction - Under Construction - Under Construction - Under Construction - Under Construction**

## I. Classifying ITS1 reads
### Preprocess ITS1 reads
This workflow assumes you have one or more fasta files of ITS1 reads amplified with 5'-GTAAAAGTCGTAACAAGGTTTC and 5'-GTTCAAAGAYTCGATGATTCAC (18SF,5.8SR1). This workflow uses mothur (http://mothur.org/) for some of the steps. If you have fastq files, you'll need to convert them to fastas (using, for instance, mothur's fastq.info command) then chop the reads to 200bp and subsample to 5000 reads. The code below assumes you are using a SLURM batch system. The make_swarm.pl command is a helper command for formatting swarm submissions. In the instructions below, it is assumed that you've set the $HFT variable to the path for this repository, for example: export HFT="../../fungal_toolkit_XH"

````
module load mothur/1.39.1
perl -w $HFT/make_swarm.pl "mothur '#make.group(fasta=FOO,groups=BAR);sub.sample(fasta=current,group=current,size=5000);chop.seqs(fasta=current,numbases=200,short=T,processors=4)'" -- *.fasta
...
swarm -f swarm_file -t 4 -g 4 --time 24:00:00 --no-comment --module mothur
...
````

### Concatenate all the subsampled reads and run blast (2.11.0+) against both the ITS1db4 and the Candida speciation database
Note that you will need to format the ITS and Candida database files before use. Use something like: makeblastdb -in file.fa -dbtype nucl
````
cat *.subsample.chop.fasta > input.ITS.sub5000.chop.fasta
cat *.subsample.groups > input.ITS.sub5000.groups
module load blast
blastn -task blastn -db $HFT/ITSdb4M.mod.fasta -query input.ITS.sub5000.chop.fasta -out ITS_sub5000.v.ITSdb4.mod.blastn.out -evalue 1e-6 -outfmt '6 std qcovs' -num_threads 32 -max_target_seqs 10
blastn -task blastn -db $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa -query input.ITS.sub5000.chop.fasta -out ITS_sub5000.v.CandidaITS1.blastn.out -evalue 1e-6 -outfmt 6 -num_threads 32 -max_target_seqs 10
````

#### Database: ITSdb4
This is a rebuilt version of the ITS database from Findley et al (https://pubmed.ncbi.nlm.nih.gov/23698366/). Full description in Huang et al.
````
09:48 biowulf fungal_toolkit_XH$ grep -c '>' ITSdb4M.mod.fasta 
54975
09:48 biowulf fungal_toolkit_XH$ grep '>' ITSdb4M.mod.fasta | head -n 3
>1001844_Root;Fungi;Basidiomycota;unclassified_Basidiomycota_class;unclassified_Basidiomycota_order;unclassified_Basidiomycota_family;unclassified_Basidiomycota_genus
>1002753_Root;Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricales_family_incertae_sedis;Lentinula
>1002759_Root;Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricales_family_incertae_sedis;Lentinula
````

#### Database: Candida_species_ITS1
This is a database of Candida species ITS1 sequences, including Candida auris clade I-IV
````
13:22 biowulf fungal_toolkit_XH$ grep -c '>' Candida_species_ITS1.trimmed.cdhit1.0.fa
89
13:23 biowulf fungal_toolkit_XH$ grep '>' Candida_species_ITS1.trimmed.cdhit1.0.fa | head -n 3
>AB040215.1_Candida_albicans_
>AB365317.1_Candida_albicans_Ethiopia
>AB369915.1_Candida_albicans_Brazil
````

### Parsing BLAST output
Genus-level classification (60% query coverage, 70% identity)
````
usage: blastn_parsing_script_ITSdb4.pl blast_output query_fasta coverage_cutoff %id_cutoff output
````
The following commands parse the ITSdb4 BLAST output and then turn that into a taxonomy table at the genus-level using mothur. tax.summary.parser.pl is a helper script that simplifies the mothur tax.summary file for plotting in Excel/R.
````
perl $HFT/blastn_parsing_script_ITSdb4.pl ITS_sub5000.v.ITSdb4.mod.blastn.out input.ITS.sub5000.chop.fasta 60 70 ITS_sub5000.v.ITSdb4.mod.blastn.taxonomy
mothur "#summary.tax(taxonomy=ITS_sub5000.v.ITSdb4.mod.blastn.taxonomy,group=input.ITS.sub5000.groups)"
perl -w $HFT/tax.summary.parser.pl -in ITS_sub5000.v.ITSdb4.mod.blastn.tax.summary -prepend_phylum -taxlevel genus -maxtax 25 > ITS_ITS1_sub5000_blastn_cov60pid70_genus.txt
````

Candida Species-level classification (50 nt HSP, 95% identity)
````
usage: blastn_parsing_script_CandidaITS1_cov90.pl blast_output query_fasta aln_length_cutoff %id_cutoff output
````
The following commands parse the Candida ITS BLAST output and then turn that into a taxonomy table at the species-level using mothur:
````
perl $HFT/blastn_parsing_script_CandidaITS1_cov90.pl ITS_sub5000.v.CandidaITS1.blastn.out input.ITS.sub5000.chop.fasta 50 95 ITS_sub5000.v.CandidaITS1.blastn.taxonomy
mothur "#summary.tax(taxonomy=ITS_sub5000.v.CandidaITS1.blastn.taxonomy,group=input.ITS.sub5000.groups)"
perl -w $HFT/tax.summary.parser.pl -in ITS_sub5000.v.CandidaITS1.blastn.tax.summary -prepend_phylum -taxlevel species -maxtax 12 > ITS_ITS1_sub5000_blastn_species.txt
````

Get C. auris clade information for the samples. Produces a table of read counts by C. auris clade
````
usage: blastn_parsing_script_CandidaITS1_cov90_strains.pl blast_output aln_length_cutoff groups_file %id_cutoff output
````
Parse BLAST output into clade-level calls (70 nt HSP, 95% identity)
````
perl $HFT/blastn_parsing_script_CandidaITS1_cov90_strains.pl ITS_sub5000.v.CandidaITS1.blastn.out 70 95 input.ITS.sub5000.groups Cauris_clades.out
````

## II. Testing the Candida speciation database

Compare the Candida speciation database to itself. Parse out inter/intra-species hits.
````
blastn -task blastn -db $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa -query $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa -out Candida_species_ITS1.trimmed.cdhit1.0_vs_self_blastn.out -evalue 1e-6 -outfmt 6 -num_threads 4
perl $HFT/blastn_parsing_script_ITS1db_test.pl Candida_species_ITS1.trimmed.cdhit1.0_vs_self_blastn.out $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa Candida_species_ITS1.trimmed.cdhit1.0_blastn_intra_parsed Candida_species_ITS1.trimmed.cdhit1.0_blastn_inter_parsed
````
Add 1% sequencing errors and rerun blast. The error rate can be changed by editing the IntroduceRandomError.pl script. Parse out inter/intra-species hits.
````
perl $HFT/IntroduceRandomError.pl $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa Candida_species_ITS1.trimmed.cdhit1.0.fa.ErrorAdded
blastn -task blastn -db $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa -query Candida_species_ITS1.trimmed.cdhit1.0.fa.ErrorAdded -out Candida_species_ITS1.trimmed.cdhit1.0_1percentErrorAdded_vs_self_blastn.out -evalue 1e-6 -outfmt 6 -num_threads 4
perl $HFT/blastn_parsing_script_ITS1db_ErrorAddedTest.pl Candida_species_ITS1.trimmed.cdhit1.0_1percentErrorAdded_vs_self_blastn.out $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa Candida_species_ITS1.trimmed.cdhit1.0_1percentErrorAdded_blastn_intra_parsed Candida_species_ITS1.trimmed.cdhit1.0_1percentErrorAdded_blastn_inter_parsed
````
Summarize the number in intra- vs inter-species hits for the top BLAST hits. There shouldn't be any hits not within the same species, but the number of hits can exceed the number of queries if there are multiple best hits with identical bit scores.
````
#1% error run
perl $HFT/blastn_parsing_script_ITS1db_ErrorAddedTest_CountBestHits.pl Candida_species_ITS1.trimmed.cdhit1.0_1percentErrorAdded_vs_self_blastn.out $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa
The number of best hits that are within the same species: 89
The number of best hits that are not within the same species: 0
#identity run (no errors)
perl $HFT/blastn_parsing_script_ITS1db_ErrorAddedTest_CountBestHits.pl Candida_species_ITS1.trimmed.cdhit1.0_vs_self_blastn.out $HFT/Candida_species_ITS1.trimmed.cdhit1.0.fa
The number of best hits that are within the same species: 89
The number of best hits that are not within the same species: 0
````
