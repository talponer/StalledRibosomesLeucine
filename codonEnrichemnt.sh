######################################################################
### This script generates genome browser tracks and calculate single
### codon occupancy:

# Add the modules:
module add samtools/1.9;
module add fastx_toolkit/0.0.14;
module add bowtie2/2.3.5;
module add STAR/2.7.0f
module add cutadapt/3.5;
module add chip-seq/1.5.5
module add bedtools/2.28.0
module add UMI-tools/1.0.0git
module add getLIMS/6.0.0
module add pipeline/0.3.2
module add rpfTools/1.2.0
module add ucsc-tools/dec-10-2019


######################################################################
### Merge sequencing runs. I will merge library PF35 and PF36 since
### they have the same samples:
for I in $(grep "PF35" 2022-05-03_samples.txt | awk 'BEGIN{FS=OFS="_"}{print $1}' | uniq )
do
    for K in $(cut -f1 umi_data/${I}_51_1_001_whitelist.pipeline)
    do
	## SID=$(echo $I | awk 'BEGIN{FS=OFS="_"}{print $1}')
	INFILE="mapping_data/sam/PF3[56]_??_?_001.filter_split_${K}.mouse_cDNA.sorted.bam"
	OUTFILE="mapping_data/sam/${I}.${K}.mouse_cDNA.sorted.bam"
	if [ ! -f $OUTFILE ]; then
	    samtools merge -@ 30 -O BAM $OUTFILE $INFILE
	fi
    done
done



######################################################################
### Make genome browser tracks

## 1. get read counts per sample
for I in $(grep "PF35" 2022-05-03_samples.txt )
do
    for K in $(cut -f1 umi_data/${I}_whitelist.pipeline)
    do
	SID=$(echo $I | awk 'BEGIN{FS=OFS="_"}{print $1}')
	INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam"
	COV=$(samtools view -c -F 260 $INFILE)
	echo -e "$SID\t$K\t$COV"
    done
done > 2022-05-06_samples.coverage

### Calculate resampling for each library:
## MINCOV=$(cut -f3 2022-05-06_samples.coverage | sort -k1,1n | head -1)
awk 'BEGIN{FS=OFS="\t"}{c[$0] = $3; l[$0] = $1; if(m[$1] == "" || m[$1] > $3){m[$1] = $3}}END{for(i in l){print i, m[l[i]] / c[i] }}' 2022-05-06_samples.coverage | sort -k1 -k2 > tmp
mv tmp 2022-05-06_samples.coverage

## 2. generate bedGraph files
GENOME=config/Mmusculus.GRCm38.100.cdna.transcriptLength
SHIFT=15
[ ! -d trackFiles ] && ( mkdir trackFiles )

for ID in $(cut -f1 2022-05-06_samples.coverage | sort -u )
do
    for K in $(grep $ID 2022-05-06_samples.coverage | cut -f2 )
    do
	INFILE="mapping_data/sam/${ID}.${K}.mouse_cDNA.sorted.bam"
	TMP="trackFiles/${ID}.${K}.mouse_cDNA.bg"
	OUTFILE="trackFiles/${ID}.${K}.mouse_cDNA.bw"
	RES=$(grep $ID 2022-05-06_samples.coverage | grep $K | cut -f4)
	echo "$ID $K $RES"
	if [ ! -f $OUTFILE ]; then
	    bamToBedGraph -i $INFILE -g $GENOME -n $RES -s $SHIFT > $TMP
	    bedGraphToBigWig $TMP $GENOME $OUTFILE
	    rm $TMP
	fi
    done
done

######################################################################
### Calculate single codon coverage
# 0. Make SGA files from BAM (predicted A-site)
[ ! -d mapping_data/sga ] && ( mkdir mapping_data/sga )

for SID in $(grep "PF35" 2022-05-06_samples.coverage | cut -f1 | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	echo "  $K"
	INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam"
	OUTFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	samtools sort -@ 25 -n $INFILE 2> /dev/null | bamToSga -i /dev/stdin -s 15 > $OUTFILE
    done
done

## 1. count reads in codons
[ ! -d codonCounts ] && ( mkdir codonCounts )

for SID in $(grep "PF35" 2022-05-06_samples.coverage | cut -f1 | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	INFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	REF="codonCounts/2022-05-06.mouse_cDNA.codonCounts.dat"
	OUTFILE="codonCounts/tempCounts.dat"
	[ ! -f $REF ] && ( cp config/Mus_musculus.GRCm38.100.codons.APPRIS_uniq.sga $REF)
	sort -k1,1 -k3,3n $REF $INFILE |
	    chipscore -c 99999999 -A "cod" -B "cDNA" -b 0 -e 2 -t 0 > $OUTFILE
	mv $OUTFILE $REF
	echo -e "$SID\t$K" >> 2022-05-06_samples.order
    done
done

## 2. remove empty codons:
awk '{c = 0; for(i=6; i<=NF; i++) {c+=$i}; if(c > 8){print}}' codonCounts/2022-05-06.mouse_cDNA.codonCounts.dat > codonCounts/2022-05-06.mouse_cDNA.codonCountsShort.dat

## 3. change counts file for downstream analyses
awk 'BEGIN{FS=OFS="\t"; while( (getline < "config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = id[1]}}{print n[$1], $1, $3, $6,$7,$8,$9,$10,$11,$12,$13}' codonCounts/2022-05-06.mouse_cDNA.codonCountsShort.dat > codonCounts/2022-05-06.mouse_cDNA.codonCountsMod.dat

### Enlarge the region to span 3 codons
for SID in $(cut -f1 2022-05-06_samples.coverage | grep "PF35" | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	INFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	REF="codonCounts/2022-05-09.mouse_cDNA.codonCounts.dat"
	OUTFILE="codonCounts/tempCounts.dat"
	[ ! -f $REF ] && ( cp config/Mus_musculus.GRCm38.100.codons.APPRIS_uniq.sga $REF)
	sort -k1,1 -k3,3n $REF $INFILE |
	    chipscore -c 99999999 -A "cod" -B "cDNA" -b -3 -e 5 -t 0 > $OUTFILE
	mv $OUTFILE $REF
    done
done

## remove empty codons:
awk '{c = 0; for(i=6; i<=NF; i++) {c+=$i}; if(c > 8){print}}' codonCounts/2022-05-09.mouse_cDNA.codonCounts.dat > codonCounts/2022-05-09.mouse_cDNA.codonCountsShort.dat

## change counts file for downstram analyses
awk 'BEGIN{FS=OFS="\t"; while( (getline < "config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = id[1]}}{print n[$1], $1, $3, $6,$7,$8,$9,$10,$11,$12,$13}' codonCounts/2022-05-09.mouse_cDNA.codonCountsShort.dat > codonCounts/2022-05-09.mouse_cDNA.codonCountsMod.dat

######################################################################
### countig RPF reads per transcript
ANNO=config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed
[ ! -d readCounts ] && ( mkdir readCounts )

for SID in $(cut -f1 2022-05-06_samples.coverage | grep "PF35" | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	echo "  $K"
	INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam"
	OUTFILE="readCounts/${SID}.${K}.readCounts.dat"
	LOGFILE="readCounts/${SID}.${K}.readCounts.log"
	samtools sort -@ 25 -n $INFILE 2> /dev/null | bamToBed -i - | rpf-counts -a $ANNO -s 15 > $OUTFILE 2> $LOGFILE
    done
done



######################################################################
### Get read enrichment near WT codns that change occupation in
### deficient diet

PEAK1=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.sga
PEAK2=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.sga
PEAK3=analysis/results/2022-05-09/ZAKdeficientVsStandardDiet_differentialCodons_UP.sga
PEAK4=analysis/results/2022-05-09/ZAKdeficientVsStandardDiet_differentialCodons_DOWN.sga

[ ! -d chipCor/2022-05-06 ] && ( mkdir chipCor/2022-05-06 )

for SID in $(grep 'PF35' 2022-05-06_samples.coverage | cut -f1 | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	echo "  $K"
	INFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	OUTFILE1="chipCor/2022-05-06/${SID}.${K}_vs_WT_DD.Enriched.dat"
	OUTFILE2="chipCor/2022-05-06/${SID}.${K}_vs_WT_DD.Depleted.dat"
	cat $INFILE $PEAK1 |
	    sort -k1,1 -k3,3n |
	    chipcor -A "peak" \
		    -B "cDNA" \
		    -b -99 \
		    -e 100 \
		    -w 1 \
		    -c 10000 \
		    -n 0 > $OUTFILE1
	cat $INFILE $PEAK2 |
	    sort -k1,1 -k3,3n |
	    chipcor -A "peak" \
		    -B "cDNA" \
		    -b -99 \
		    -e 100 \
		    -w 1 \
		    -c 10000 \
		    -n 0 > $OUTFILE2
    done
done


### Use chip-extract to do the same
for SID in $(grep 'PF35' 2022-05-06_samples.coverage | cut -f1 | sort -u )
do
    echo $SID
    for K in $(grep $SID 2022-05-06_samples.coverage | cut -f2 )
    do
	echo "  $K"
	INFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	OUTFILE1="chipCor/2022-05-06/${SID}.${K}_vs_WT_DD.Enriched.mat"
	OUTFILE2="chipCor/2022-05-06/${SID}.${K}_vs_WT_DD.Depleted.mat"
	cat $INFILE $PEAK1 |
	    sort -k1,1 -k3,3n |
	    chipextract -A "peak" \
			-B "cDNA" \
			-b -99 \
			-e 100 \
			-w 1 \
			-c 10000 \
			> $OUTFILE1
	cat $INFILE $PEAK2 |
	    sort -k1,1 -k3,3n |
	    chipextract -A "peak" \
			-B "cDNA" \
			-b -99 \
			-e 100 \
			-w 1 \
			-c 10000 \
			> $OUTFILE2
    done
done

## check distribution of up and down codons:
sed 's/peak/pUp/' $PEAK1 |
    sort -k1,1 -k3,3n $PEAK2 - |
    chipcor -A "pUp" \
	    -B "peak" \
	    -b -9999 \
	    -e 10000 \
	    -w 50 \
	    -c 10000 \
	    -n 0 > chipCor/2022-05-06/WT_DD_UP_vs_DOWN.dat


## Now for the mutant:
sed 's/peak/pUp/' $PEAK3 |
    sort -k1,1 -k3,3n $PEAK4 - |
    chipcor -A "pUp" \
	    -B "peak" \
	    -b -9999 \
	    -e 10000 \
	    -w 50 \
	    -c 10000 \
	    -n 0 > chipCor/2022-05-06/ZAK_DD_UP_vs_DOWN.dat



# sed 's/peak/pUp/' $PEAK3 |
#     sort -k1,1 -k3,3n $PEAK4 - |
#     chipcor -A "pUp" \
# 	    -B "peak" \
# 	    -b -9999 \
# 	    -e 10000 \
# 	    -w 50 \
# 	    -c 10000 \
# 	    -n 0 > chipCor/HFD.NAC_Chow.NAC_UP_vs_DOWN.dat

## Autocorrelation plot for UP and DOWN regions:
chipcor -A "peak" \
	-B "peak" \
	-b -2999 \
	-e 3000 \
	-w 3 \
	-c 10000 \
	-n 0 $PEAK1 > chipCor/2022-05-06/WT_DD_UP_autocorrelation.dat

chipcor -A "peak" \
	-B "peak" \
	-b -2999 \
	-e 3000 \
	-w 3 \
	-c 10000 \
	-n 0 $PEAK2 > chipCor/2022-05-06/WT_DD_DOWN_autocorrelation.dat

# chipcor -A "peak" \
# 	-B "peak" \
# 	-b -2999 \
# 	-e 3000 \
# 	-w 3 \
# 	-c 10000 \
# 	-n 0 $PEAK3 > chipCor/HFD.NAC_UP_autocorrelation.dat

######################################################################
### Find AA under the pausing sites

## 1. make bed files:
## awk 'BEGIN{FS=OFS="\t"; while( (getline < "config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = $1}}$3 > 31{print n[$1], $3-31, $3+32} ' $PEAK1 | awk 'BEGIN{c = 0; g = ""; FS=OFS="\t"} {if (g == ""){g = $1}; if($2 > c && g == $1){print}; if ($1 == g) {c = $2+3}else{print; c = $2+3; g = $1}}' > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.bed
awk 'BEGIN{FS=OFS="\t"; while( (getline < "config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = $1}}$3 > 31{print n[$1], $3-31, $3+32} ' $PEAK1 > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.bed

## 2. get the FASTA from the data
BED1=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.bed
FASTA=/data/databases/mouse/fasta/Mmusculus.GRCm38.100.cdna.ensembl.fa

bedtools getfasta -tab -fi $FASTA -bed $BED1 |
    fastaToAA.py |
    sed '/Ter/d' |
    sort -u > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP_AminoA.dat

## 3. make a pfm:
cat analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP_AminoA.dat | makePwm.py > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP_AminoA.pwm

# cat analysis/results/2022-05-09/HFD.NAC_Chow.NAC_differentialCodons_UP_AminoA.dat | makePwm.py > analysis/results/2022-05-09/HFD.NAC_Chow.NAC_differentialCodons_UP_AminoA.pwm


### Now for the down:
awk 'BEGIN{FS=OFS="\t"; while( (getline < "config/Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = $1}}$3 > 31{print n[$1], $3-31, $3+32} ' $PEAK2 > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.bed

## 2. get the FASTA from the data
BED2=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.bed

bedtools getfasta -tab -fi $FASTA -bed $BED2 |
    fastaToAA.py |
    sed '/Ter/d' |
    sort -u > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN_AminoA.dat


## calculate total number of AA in the whole region
bedtools getfasta -tab -fi $FASTA -bed $BED1 |
    fastaToAA.py |
    sed '/Ter/d' |
    uniq |
    awk '$1 !~ ">"{split($1, n, ""); for(i=8; i<=length(n)-7; i++){c[n[i]]++}}END{for(i in c){print i, c[i]}}' |
    sort -k2,2n > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.aaFreq

bedtools getfasta -tab -fi $FASTA -bed $BED2 |
    fastaToAA.py |
    sed '/Ter/d' |
    uniq |
    awk '$1 !~ ">"{split($1, n, ""); for(i=8; i<=length(n)-7; i++){c[n[i]]++}}END{for(i in c){print i, c[i]}}' |
    sort -k2,2n > analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.aaFreq

## Calculate BG aa frequency
BED3=config/Mus_musculus.GRCm38.100.startStop.APPRIS_expressed.bed
bedtools getfasta -tab -fi $FASTA -bed $BED3 |
    fastaToAA.py |
    sed '/Ter/d' |
    uniq |
    awk '$1 !~ ">"{split($1, n, ""); for(i=1; i<=length(n); i++){c[n[i]]++}}END{for(i in c){print i, c[i]}}' |
    sort -k2,2n > config/Mus_musculus.GRCm38.100.startStop.APPRIS_expressed.aaFreq

######################################################################
### Get the paks that have a Y in the center (+-3 aa) and make an
### auto-correlation plot to see if disome signature is there
######################################################################
PEAKY1=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP_Y.sga
PEAKY2=analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN_Y.sga

bedtools getfasta -tab -fi $FASTA -bed $BED1 |
    fastaToAA.py |
    awk '{print substr($1, 8, 7)}' |
    awk '$1 ~ "Y" {print NR}' |
    awk 'BEGIN{FS=OFS="\t"; c=1; while( (getline < "analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_UP.sga") >0 ){sga[c] = $0; c++}}{print sga[$1]}' > $PEAKY1

bedtools getfasta -tab -fi $FASTA -bed $BED2 |
    fastaToAA.py |
    awk '{print substr($1, 8, 7)}' |
    awk '$1 ~ "Y" {print NR}' |
    awk 'BEGIN{FS=OFS="\t"; c=1; while( (getline < "analysis/results/2022-05-09/WTdeficientVsStandardDiet_differentialCodons_DOWN.sga") >0 ){sga[c] = $0; c++}}{print sga[$1]}' > $PEAKY2


chipcor -A "peak" \
	-B "peak" \
	-b -2999 \
	-e 3000 \
	-w 3 \
	-c 10000 \
	-n 0 $PEAKY1 > chipCor/2022-05-06/WT_DD_UP_Y_autocorrelation.dat

chipcor -A "peak" \
	-B "peak" \
	-b -2999 \
	-e 3000 \
	-w 3 \
	-c 10000 \
	-n 0 $PEAKY2 > chipCor/2022-05-06/WT_DD_DOWN_Y_autocorrelation.dat


sed 's/peak/pUp/' $PEAK1 |
    sort -k1,1 -k3,3n $PEAKY1 - |
    chipcor -A "peak" \
	    -B "pUp" \
	    -b -2999 \
	    -e 3000 \
	    -w 3 \
	    -c 1 \
	    -n 0 > chipCor/2022-05-06/WT_DD_UP_Y_corr.dat

