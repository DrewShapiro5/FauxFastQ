python3 ~/Desktop/Lab/FauxFastQ/FauxFQ.py sperm.referenceseq.fa o_f.fq o_r.fq -n 1 --reference_forward reference_reads/read1.fastq.gz --reference_reverse reference_reads/read2.fastq.gz -d -i -r --append
xattr -d com.apple.quarantine o_f.fq
xattr -d com.apple.quarantine o_r.fq
