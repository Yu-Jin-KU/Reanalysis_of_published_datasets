##suppa & salmon

suppa.py generateEvents -i atRTD3_TS_21Feb22_transfix.gtf  -o ATRTD3.events -e SE SS MX RI FL -f ioe
$ wc -l *ioe


awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' *.ioe > ATRTD3.all.events.ioe
wc -l ATRTD3.all.events.ioe

salmon index -t atRTD3_29122021.fa -i transcripts_index -k 31

##salmon for grp7

for f1 in *_trimmed.fq.gz; do salmon quant -i transcripts_index -p 8 --libType A   -r $f1 -o ${f1/_trimmed.fq.gz//}.quant;done

multipleFieldSelection.py \
-i  SRR*/quant.sf -k 1 -f 4 \
-o iso_tpm.txt

suppa.py psiPerEvent -i ATRTD3.all.events.ioe -e iso_tpm.txt -o as_events


cut -f 1-4 as_events.psi > WT_LL24.psi
cut -f 1-4 iso_tpm.txt > WT_LL24.tpm
cut -f 1,5-7 as_events.psi > WT_LL36.psi
cut -f 1,5-7 iso_tpm.txt > WT_LL36.tpm
cut -f 1,8-10 as_events.psi > GRP7_LL24.psi
cut -f 1,8-10 iso_tpm.txt > GRP7_LL24.tpm
cut -f 1,11-13 iso_tpm.txt > GRP7_LL36.tpm
cut -f 1,11-13 as_events.psi > GRP7_LL36.psi
cut -f 1,14-16 iso_tpm.txt > grp7_LL24.tpm
cut -f 1,14-16 as_events.psi > grp7_LL24.psi
cut -f 1,17-19 as_events.psi > grp7_LL36.psi
cut -f 1,17-19 iso_tpm.txt > grp7_LL36.tpm


suppa.py   diffSplice \
-m empirical -gc -i  ATRTD3.all.events.ioe   \
--save_tpm_events \
-p grp7_LL24.psi WT_LL24.psi  \
-e grp7_LL24.tpm WT_LL24.tpm \
-o grp7_LL24_diffSplice 

cat grp7_LL24_diffSplice.dpsi|perl -alne '{print if $F[2] <0.01}' |wc


suppa.py   diffSplice \
-m empirical -gc -i  ATRTD3.all.events.ioe  \
--save_tpm_events \
-p grp7_LL36.psi WT_LL36.psi  \
-e grp7_LL36.tpm WT_LL36.tpm \
-o grp7_LL36_diffSplice # 前缀

cat grp7_LL36_diffSplice.dpsi|perl -alne '{print if $F[2] <0.01}' |wc


