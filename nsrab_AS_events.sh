##suppa & salmon

suppa.py generateEvents -i atRTD3_TS_21Feb22_transfix.gtf  -o ATRTD3.events -e SE SS MX RI FL -f ioe
$ wc -l *ioe


awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > ATRTD3.all.events.ioe
wc -l ATRTD3.all.events.ioe

salmon index -t atRTD3_29122021.fa -i transcripts_index -k 31

##salmon for nsrab
for f1 in *_R1.fq.gz; do f2=${f1/_R1/_R2} && echo $f1 $f2 && salmon quant -i transcripts_index -p 6 --libType A -1 $f1 -2 $f2 -o ${f1/_R1.fastq.gz//}.quant;done

multipleFieldSelection.py \
-i  SRR*/quant.sf -k 1 -f 4 \
-o iso_tpm.txt

suppa.py psiPerEvent -i ATRTD3.all.events.ioe -e iso_tpm.txt -o as_events


cut -f 1,8-10 as_events.psi > WT.psi
cut -f 1,8-10 iso_tpm.txt > WT.tpm
cut -f 1,11-13 as_events.psi > nsrab.psi
cut -f 1,11-13 iso_tpm.txt > nsrab.tpm




suppa.py   diffSplice \
-m empirical -gc -i  ATRTD3.all.events.ioe   \
--save_tpm_events \
-p nsrab.psi WT.psi  \
-e nsrab.tpm WT.tpm \
-o nsrab_diffSplice 

cat project_diffSplice.dpsi|perl -alne '{print if $F[2] <0.01}' |wc
