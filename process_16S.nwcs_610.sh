### nwcs610+ZEBS BMK project (16S breast milk sequenced at CHOP Microbiome Facility) ###

############################################################################
### 1. Start from per-sample FASTQ files for DADA2 and phyloseq analysis
### a. Run DADA2 and BLAST for species assignment
#	INPUT:
#		220525_M03543_0295_000000000-KDCGT_L001
#		220525_M04734_0349_000000000-KDCGW_L001
#		220527_M03543_0296_000000000-KCRV2_L001
#		220531_M03543_0297_000000000-KCR4P_L001
#		nwcs610_ZEBS_combined_batching.samples.052421.txt
#	EXECUTION:
OUTDIR=/Lab_Share/PROMISE/nwcs610
mkdir -p ${OUTDIR}/phyloseq

# remove empty gzipped fastq files
find . -type f -name "*.fastq.gz" -size -50c -exec rm -i {} \;
find . -type f -name "*.fastq.gz" -size -50c -exec rm -i {} \;

./run_DADA2.nwcs_610.R

/share/BLAST/bin/blastn -task blastn -query ${OUTDIR}/fastq/DADA2_sequences.fasta -db /Lab_Share/SILVA/SILVA_138_SSURef_NR99_tax_silva_trunc -out ${OUTDIR}/fastq/BLAST_results.txt -num_threads 24 -evalue 0.01 -outfmt 6 -max_target_seqs 10
ruby parse_blast_output.rb --type SILVA --ignore_strain --consensus weighted ${OUTDIR}/fastq/BLAST_results.txt /Lab_Share/SILVA/SILVA_138_SSURef_NR99_tax_silva_trunc.fixed_headers > ${OUTDIR}/fastq/BLAST_results.parsed.txt

### b. phyloseq and such


#### c. Trying LEfSe
#lefse-format_input.py ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.txt ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.in -c 1 -o 1000000
#run_lefse.py ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.in ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.res
#lefse-plot_res.py --format pdf ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.res ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.pdf

#lefse-format_input.py ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.txt ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.in -c 1 -o 1000000
#run_lefse.py ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.in ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.res
#lefse-plot_res.py --format pdf ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.res ${OUTDIR}/phyloseq/CHANGES_decontam_filtered_for_LEfSe.Cohort.pdf



