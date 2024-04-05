
make_minigenome:
	../make_mini_genome.pl  --gtf data/annotation.gtf --genome_fa data/genome.fa --shrink_introns --accs_restrict_file genes.want.txt


make_minigenome_from_transcript_list:
	../make_mini_genome.pl  --gtf data/annotation.gtf --genome_fa data/genome.fa --shrink_introns --accs_restrict_file transcripts.want.txt --out_prefix mini_by_transcripts


trim_utrs: make_minigenome
	../util/gtf_UTR_trimmer.pl  minigenome.gtf minigenome.fa  > minigenome.UTRs_trimmed.gtf

extract_minigenome_align_reads:
	hisat2_extract_splice_sites.py minigenome.gtf  > minigenome.gtf.ss
	hisat2_extract_exons.py  minigenome.gtf > minigenome.gtf.exons
	hisat2-build --exon ./minigenome.gtf.exons --ss ./minigenome.gtf.ss  ./minigenome.fa ./minigenome.fa
	hisat2 --dta -x minigenome.fa -1 data/reads_1.fq -2 data/reads_2.fq | samtools view -Sb -F 4 | samtools sort -o alignments.hisat2.bam
	../util/get_aligned_reads_fastqs.pl  --bam alignments.hisat2.bam --left_fq data/reads_1.fq --right_fq data/reads_2.fq --output_prefix mini

clean:
	rm -f ./mini* alignments.hisat2.bam

