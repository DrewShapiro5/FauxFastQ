python3 ../FauxFQ.py sperm.referenceseq.fa o5_f.fq o5_r.fq -n 5 --reference_forward reference_reads/read1.fastq.gz --reference_reverse reference_reads/read2.fastq.gz -d -i -r
xattr -d com.apple.quarantine o_f.fq
xattr -d com.apple.quarantine o_r.fq
