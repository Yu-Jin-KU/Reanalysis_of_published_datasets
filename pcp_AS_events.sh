##suppa & salmon

suppa.py generateEvents -i atRTD3_TS_21Feb22_transfix.gtf  -o ATRTD3.events -e SE SS MX RI FL -f ioe
$ wc -l *ioe


awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > ATRTD3.all.events.ioe
wc -l ATRTD3.all.events.ioe

salmon index -t atRTD3_29122021.fa -i transcripts_index -k 31

##salmon for pcp
for f1 in *_paired_R1.fastq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && salmon quant -i transcripts_index/ -p 6 --libType A -1 $f1 -2 $f2 -o ${f1/_R1.fastq.gz//}.quant;done

multipleFieldSelection.py -i  ../*_paired*/.quant/quant.sf -k 1 -f 4 -o iso_tpm.txt

suppa.py psiPerEvent -i ATRTD3.all.events.ioe -e iso_tpm.txt -o as_events

cut -f 1-2 as_events.psi > WT.psi
cut -f 1-2 iso_tpm > WT.tpm
cut -f 1,3 as_events.psi > pcp.psi
cut -f 1,3 iso_tpm > pcp.tpm


suppa.py   diffSplice  -m empirical -gc -i  ATRTD3.all.events.ioe --save_tpm_events  -p pcp.psi WT.psi  -e pcp.tpm WT.tpm  -o pcp_diffSplice
