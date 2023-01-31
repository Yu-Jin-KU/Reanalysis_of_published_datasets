##suppa & salmon

suppa.py generateEvents -i atRTD3_TS_21Feb22_transfix.gtf  -o ATRTD3.events -e SE SS MX RI FL -f ioe
$ wc -l *ioe


awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > ATRTD3.all.events.ioe
wc -l ATRTD3.all.events.ioe

salmon index -t atRTD3_29122021.fa -i transcripts_index -k 31


##salmon for flail
for f1 in *_paired_R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && salmon quant -i transcripts_index -p 6 --libType A -1 $f1 -2 $f2 -o ${f1/_R1.fastq.gz//}.quant;done
multipleFieldSelection.py -i  RTD*/.quant/quant.sf -k 1 -f 4 -o iso_tpm.txt

suppa.py psiPerEvent -i ATRTD3.all.events.ioe -e iso_tpm.txt -o as_events2

cut -f 1-4 as_events2.psi > flail.psi
cut -f 1-4 iso_tpm.txt > flail.tpm
cut -f 1,5-7 as_events2.psi > WT.psi
cut -f 1,5-7 iso_tpm.txt > WT.tpm


suppa.py   diffSplice \
-m empirical -gc -i  ATRTD3.all.events.ioe  \
--save_tpm_events \
-p flail.psi WT.psi   \
-e flail.tpm WT.tpm \
-o flail_diffSplice   
