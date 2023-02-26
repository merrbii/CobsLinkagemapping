# Pipeline used for linkage mapping analysis in _Cardiocondyla obscurior_
##### By Mohammed Errbii (Simo)

Here sampled and sequenced a total of 104 individuals from a cross between a queen and a male of the ant _Cardiocondyla obscurior_. This included two grandparents, two F1 parents and 100 F2 individuals (88 workers and 12 males). Details about the cross, our questions as well as our main findings can be found [here]().
The genome assembly and annotation of _C. obscurior_ as well as the raw reads can be downloaded from NCBI under BioProject: PRJNA934066

### Filter raw reads using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```bash
ls *_R1_001.fastq.gz|sed 's/_R1_001.fastq.gz//g'|nice parallel --jobs 15 --eta 'java -jar /global/projects/programs/bin/trimmomatic-0.38.jar PE {}_R1_001.fastq.gz {}_R2_001.fastq.gz paired/{}_paired_1.fastq.gz unpaired/{}_unpaired_1.fastq.gz paired/{}_paired_2.fastq.gz unpaired/{}_unpaired_2.fastq.gz SLIDINGWINDOW:4:20 MINLEN:40'
```

### Map trimmed reads with [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)

```bash

for infile in *_paired_1.fastq.gz; do base=$(basename -s _paired_1.fastq.gz ${infile}); bwa mem -t 25 ~/genome/Cobs3.1.clean.fa ${infile} ${base}_paired_2.fastq.gz > samfiles/${base}.sam;done
```

### Pre-processing alignments using [Picard](https://broadinstitute.github.io/picard/) and [SAMtools](http://samtools.sourceforge.net/)

```bash
# a- clean sam files
ls *sam|parallel --progress --eta --jobs 5 'java -jar /global/projects/programs/bin/picard.jar CleanSam I="{}" O="{=s/(.*)\..*/$1/=}C.sam"'

# b- covert to bam and sort bam file (sort by default assume sorting by coordinates):
ls *C.sam| nice parallel --progress --eta --jobs 5  'samtools view -S -u -b {} | samtools sort -o {.}S.bam'

# c- fix mate information if needed
ls *CS.bam| nice parallel --progress --eta --jobs 5  'java -jar /global/projects/programs/bin/picard.jar FixMateInformation I="{}" O="{.}F.bam" ADD_MATE_CIGAR=true ASSUME_SORTED=true'

# d- fix read group information
nice parallel -a rgid.txt --colsep '\t' --jobs 20 --eta java -jar /global/projects/programs/bin/picard.jar AddOrReplaceReadGroups I={1} O={5}R.bam RGID={2} RGLB={4} RGPL=illumina RGPU={3} RGSM={4}
# assumes you have a tab separated rgid file with following columns: <bamfile name> <readgroup id>  <lane if any> <sample name (unique)>  <basename for output>

# e- locates and tags duplicate reads in bam file
ls *R.bam |nice parallel --progress --eta --jobs 5 'java -jar /global/projects/programs/bin/picard.jar MarkDuplicates I="{}" O="{.}D.bam" M="{.}D.txt" MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ASSUME_SORTED=true'

# f- index bam file
samtools index alignment.CSFRD.bam
```


### Variant calling using [GATK](https://gatk.broadinstitute.org/hc/en-us)

```bash

# a- generate gVCFs
ls *D.bam |parallel --progress --eta --jobs 10 '/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar /global/projects/programs/source/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar HaplotypeCaller -R ~/genome/Cobs3.1.clean.fa -I {} -O {.}.g.vcf.gz -ERC GVCF'

# b- combine gVCFs
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar /global/projects/programs/source/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar CombineGVCFs -R ~/genome/Cobs3.1.clean.fa -V input.list -O Cobs3.1.linkageMapping.g.vcf.gz

# c- genotype the combined gVCF
/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -jar /global/projects/programs/source/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar GenotypeGVCFs -R ~/genome/Cobs3.1.clean.fa -V Cobs3.1.linkageMapping.g.vcf.gz -O Cobs3.1.linkageMapping.vcf.gz
```

### Variant filtering


Get bi allelic SNPs
```bash

vcftools --gzvcf gCobs3.1.linkageMapping.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out Cobs3.1.linkageMapping.snps

# compress and index
bgzip Cobs3.1.linkageMapping.snps.recode.vcf
tabix Cobs3.1.linkageMapping.snps.recode.vcf.gz
```
Get heterozygous markers in F1 queen (F1Queen) that are homozygous fixed for different alleles in the parents or heterozygous in grandmother (PQueen)
```bash
gatk SelectVariants -V Cobs3.1.linkageMapping.snps.recode.vcf.gz -O biall.snps.noMendelianErrParents.vcf.gz -select 'vc.getGenotype("F1Queen").isHet() && (vc.getGenotype("PQueen").isHomRef() && vc.getGenotype("PMale").isHomVar() || vc.getGenotype("PQueen").isHet() && vc.getGenotype("PMale").isHom()) || vc.getGenotype("F1Queen").isHet() && (vc.getGenotype("PQueen").isHomVar() && vc.getGenotype("PMale").isHomRef() || vc.getGenotype("PQueen").isHet() && vc.getGenotype("PMale").isHom())' -select 'vc.getGenotype("F1Male").isHomVar() && (vc.getGenotype("PQueen").isHomVar() || vc.getGenotype("PQueen").isHet()) || vc.getGenotype("F1Male").isHomRef() && (vc.getGenotype("PQueen").isHomRef() || vc.getGenotype("PQueen").isHet())'
```
In ants males are haploid. Any position called as heterozygous in males indicates a genotyping error and should be excluded. get only positions where none of the males is heterozygous
```bash
gatk SelectVariants -V biall.snps.noMendelianErrParents.vcf.gz -O biall.snps.noMendelianErrParents.noHetMale.vcf.gz --invertSelect true -select 'vc.getGenotype("EM1").isHet()' -select 'vc.getGenotype("EM112").isHet()'  -select 'vc.getGenotype("EM2").isHet()' -select 'vc.getGenotype("EM3").isHet()' -select 'vc.getGenotype("EM35").isHet()' -select 'vc.getGenotype("EM36").isHet()'  -select 'vc.getGenotype("EM37").isHet()' -select 'vc.getGenotype("EM38").isHet()' -select 'vc.getGenotype("EM39").isHet()' -select 'vc.getGenotype("EM4").isHet()' -select 'vc.getGenotype("EM5").isHet()' -select 'vc.getGenotype("EM78").isHet()' -select 'vc.getGenotype("F1Male").isHet()' -select 'vc.getGenotype("PMale").isHet()'
```
Check and exclude genotypes showing mendelian discrepancies using [bcftools](https://samtools.github.io/bcftools/bcftools.html) _+mendelian_ function. This concerns the diploid workers. Males get one of the two F1 queen haplotypes and we have removed heterozygous positions anyways. For the workers however, there might be instances of mendelian errors. These are either genotyping errors or novel mutations, that should be excluded. For instance if the F1 queen is 0/1 and the F1 father is 1/1 the only two possible F2 genotypes are 0/1 or 1/1 (possible genotype), but never 0/0 (**impossible genotype**). If 0/0 then that should be switched to missing. In the next section we extract the genotypes of each of the workers and compare it to the F1 father. If **impossible genotype**, we replace it with missing (./.).


```bash
## 1) get only genotypes and switch | to /
cat <(zgrep "#" biall.snps.noMendelianErrParents.noHetMale.vcf.gz) <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT[\t%GT]\n' biall.snps.noMendelianErrParents.noHetMale.vcf.gz)|sed 's@|@/@g' > biall.snps.noMendelianErrParents.noHetMale.GT.vcf

bgzip biall.snps.noMendelianErrParents.noHetMale.GT.vcf
tabix biall.snps.noMendelianErrParents.noHetMale.GT.vcf.gz

## 2)
# get worker + F1 father samples
sample=$(zgrep "^#CH" biall.snps.noMendelianErrParents.noHetMale.GT.vcf.gz|cut -f10-|tr "\t" "\n"|grep "^W")

for i in $sample; do bcftools view -s F1Male,$i biall.snps.noMendelianErrParents.noHetMale.GT.vcf.gz -Ov -o F1Male.$i.biall.snps.noMendelianErrParents.noHetMale.GT.vcf;done

## 3) loop over the different vcf files, check of genotypes of a given worker match the expected one, knowing the genotype of the father. If not replace with ./.
### replace XXX with the value from grep "^#" F1Male.$i.biall.snps.noMendelianErrParents.noHetMale.GT.vcf
for i in F1Male.W*.vcf; do base=$(basename -s .vcf $i); awk -vFS="\t" -vOFS="\t" 'NR<=XXX{print $0}; FNR > XXX {if ($10=="1/1" && $11 == "0/0") {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, "./."} else if ($10=="0/0" && $11 == "1/1") {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, "./."} else {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11} }' $i > "${base}.ed.vcf";done

## 4) get the edited worker genotypes
for i in $sample; do bcftools view -s $i F1Male.$i.biall.snps.noMendelianErrParents.noHetMale.GT.ed.vcf -Ov -o $i.biall.snps.noMendelianErrParents.noHetMale.GT.ed.vcf; done

for i in W*.vcf; do bgzip $i;done
for i in W*.vcf.gz; do tabix $i;done

## 5) get males and parents
bcftools view -s PMale,PQueen,F1Male,F1Queen,EM1,EM112,EM2,EM3,EM35,EM36,EM37,EM38,EM39,EM4,EM5,EM78 biall.snps.noMendelianErrParents.noHetMale.GT.vcf.gz -Ov -o grandParentsF1Parents.males.vcf

bgzip grandParentsF1Parents.males.vcf
tabix grandParentsF1Parents.males.vcf.gz

## 6) merge everything
bcftools merge grandParentsF1Parents.males.vcf.gz $(ls W*.ed.vcf.gz)  -Oz -o biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.vcf.gz

tabix biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.vcf.gz

## 7) clean
rm F1Male.W*
rm W*
rm grandParentsF1Parents.males.vcf.gz*

## add original annotations (i.e. metrics associated with each genotype) to the edited vcf file
bcftools annotate -a ../biall.snps.noMendelianErrParents.noHetMale.vcf.gz -c INFO,^FORMAT/GT biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.vcf.gz -Oz -o biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.vcf.gz

## 8) compare mendelian error before/after to check if everything worked
paste <(bcftools +mendelian biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.vcf.gz -T pedigree.txt ) <(bcftools +mendelian biall.snps.noMendelianErrParents.noHetMale.vcf.gz -T pedigree.txt) > afterANDbefore.txt

```
Use [VCFtools](https://vcftools.github.io/index.html) to filter based on coverage, missingness and maf.
```bash
vcftools --gzvcf biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.vcf.gz --max-meanDP 50 --min-meanDP 19 --minDP 5 --max-missing 0.8 --maf 0.2 --recode --recode-INFO-all --out biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered
```

### Phasing
Males are haploid and the direct product of the F1 queen. The 88 workers included in our design are diploid and have one haplotype from the F1 father plus the one from their F1 mother. We are looking at recombination in the F1 queen and thus interested in the maternal contribution in these workers. In the following section we exclude the father's genotype/contribution from each of the workers to only keep the maternal recombinant haplotype.

```bash
# get genotypes + position etc into a text file
cat <(bcftools view biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.vcf|grep "#C"|cut -f1-5,10-) <(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.vcf) > phasing/unphasedMarkers.txt

# replace missing (./.) genotypes with -,  homozygous alt by 1, heterozygous genotypes by 01 and homozygous ref by 0. This is what each of the sed is doing
paste <(cut -f1-5 unphasedMarkers.txt) <(cut -f6- unphasedMarkers.txt|sed 's@\./\.@-@g'|sed -E 's@1/1|1\|1@1@g'|sed -E 's@0/1|0\|1@01@g'| sed -E 's@0/0|0\|0@0@g') > unphasedMarkers.forPhasing.txt

# get F1male, F1Queen and workers' genotypes only
cut -f8,9,22- unphasedMarkers.forPhasing.txt > F1parentsF2workers.tsv

# delete parental contribution from heterozygous genotypes (01). Add these "phased" worker genotypes to the F2 males data in cut -f1-5,10-21 unphasedMarkers.forPhasing.txt and save everything into phasedMarkers.noMissF1male.tsv. exclude markers where we don't have genotype info for F1Male by awk '{if ($18 != "-") print $0}'.
paste <(cut -f1-5,10-21 unphasedMarkers.forPhasing.txt) <(awk 'BEGIN{FS=OFS="\t"}  {for (i=2; i<=NF; i++) if ($i == "01") gsub($1, "", $i)}1'  F1parentsF2workers.tsv) | awk '{if ($18 != "-") print $0}' > phasedMarkers.noMissF1male.tsv

# get annotations from unphased vcf file
bcftools view biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.vcf |grep -v "##"|awk 'BEGIN{FS=OFS="\t"}  {for (i=1; i<9; i++) printf $i FS; print "GT"}' > vcf.annotations.tsv
# change first GT to FORMAT!
```
Now switch to R to generate a phased vcf file
```R
# add annotations from vcf file before creating a vcf file
phased  <- read_delim("phasing/phasedMarkers.noMissF1male.tsv", col_names = T, delim = "\t")
ann <- read_delim("phasing/vcf.annotations.tsv", col_names = T, delim = "\t")

combined <- right_join(dplyr::select(ann,c(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)),
                       dplyr::select(phased,c(`#CHROM`,POS,ID,REF,ALT,F1Male,EM1,EM112,EM2,EM3,EM35,EM36,EM37,EM38,
                                              EM39,EM4,EM5,EM78,W101,W103,W106,W107,W109,W110,
                                              W111,W13,W14,W15,W16,W17,W18,W19,W200,W201,W202,W204,W205,W206,
                                             W207,W211,W212,W221,W23,W25,W26,W27,W29,W30,W300,W301,W302,W303,
                                             W304,W32,W33,W34,W40,W41,W42,W43,W44,W45,W46,W49,W50,W51,W52,W53,
                                             W55,W56,W57,W58,W59,W6,W60,W61,W62,W63,W64,W68,W69,W7,W70,W71,W72,
                                             W73,W74,W75,W76,W77,W78,W8,W80,W81,W82,W83,W84,W85,W86,W88,W89,W9,
                                             W90,W91,W92,W97,F1Queen)))

write.table(combined, "phasing/vcf.body.phased.tsv",
            sep = "\t",col.names = T, quote = F, row.names = F, na = "-")

```
Return to bash and produce the final clean vcf file with phased gentotypes
```bash
# add headers, "diploidize" and switch missing (-) to ./.
cat <(bcftools view biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.vcf |grep "##")  <(paste <(cut -f1-9 vcf.body.phased.tsv ) <(cat <(cut -f10- vcf.body.phased.tsv  |grep "^F") <(cut -f10- vcf.body.phased.tsv | grep -v "^F"|sed 's@1@1/1@g'|sed 's@0@0/0@g'|sed 's@-@./.@g'))) > biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.phased.vcf

bgzip biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.phased.vcf
tabix biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.phased.vcf.gz

# add annotations (in FORMAT)
bcftools annotate -a biall.snps.noMendelianErrParents.noHetMale.vcf.gz -c INFO,^FORMAT/GT phasing/biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.ANN.filtered.recode.phased.vcf.gz -Oz -o  biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.filtered.phased.vcf.gz
```
### Running [MSTmap](http://alumni.cs.ucr.edu/~yonghui/mstmap.html)


```bash
# get genotypes of the 100 F2 individuals and recode them into a format required by MSTmap, i.e. homozygous alt (1/1 into A), homozygous ref (0/0 into B) and ./. into -
bcftools view -s ^F1Male,F1Queen biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.filtered.phased.vcf.gz |bcftools query -f '%CHROM\t%POS[\t%GT]\n' |sed 's@1/1@A@g'|sed 's@0/0@B@g'|sed 's@./.@-@g'|sed 's/\./-/g' > markers.tsv

# add marker ID and sample names
cat <(bcftools view -s ^F1Male,F1Queen biall.snps.noMendelianErrParents.noHetMale.GT.noMenErr.filtered.phased.vcf.gz |grep  "#C"|cut -f1,2,10-) markers.tsv |awk 'BEGIN{FS=OFS="\t"}{$1=$1"_"$2;$2=""; print $0}'|tr -s "\t" > markers.set.for.MSTmap.txt

# then duplicate markers
awk -vOFS="\t" '{$1 = $1"-1"; print}' markers.set.for.MSTmap.txt > markers.for.MSTmap.1.tmp.txt
awk -vOFS="\t" '{$1 = $1"-2"; print}' markers.set.for.MSTmap.txt > markers.for.MSTmap.2.tmp.txt

# flip phase in one file
cat markers.for.MSTmap.2.tmp.txt|sed 's/A/0/g'|sed 's/B/A/g'|sed 's/0/B/g' > markers.for.MSTmap.2.prime.txt

# merge the two sets of markers and exclude markers on unplaced scaffolds (had less than two markers anyways)
cat markers.for.MSTmap.1.tmp.txt <(tail +2 markers.for.MSTmap.2.prime.txt)|grep -v "scaff" > markers.for.MSTmap.readyToMap.txt

# run MSTmap
~/softwares/MSTmap/mstmap markers.for.MSTmap.readyToMap.txt markers.for.MSTmap.mapped.txt >> mstmap.stdout 2>&1
```

### Manual anchoring of scaffolds
To  split chimeric scaffolds (i.e. scaffolds where markers are assigned to more than one LG), we  checked for the existence of gaps separating the parts belonging to different LGs. When gaps were absent, chimeric scaffolds were split after the last marker of the first part. For misoriented scaffolds, we visually inspected the correlation between the physical and genetic positions of markers.
```bash
# get gaps in ref:
perl -ne 'chomp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++; if($_ eq "N" && $s ==0 ){$z=$i; print "$head\t$z"; $s =1}elsif($s==1 && $_ ne "N"){$j=$i-1;print "\t$j\n";$s=0}}' ~/genome/Cobs.alpha.2.1.fa > Cobs2.1.gaps.bed

# look at composition of each LG (from MSTmap) and create a bed file with the name of scaffolds, start and end positions, and the orientation of scaffolds forming the LG.
# See attached tables with the paper for an example

# get the refined sequence for each LG using the new coordinates bed file
for i in LG*.bed; do base=$(echo $i|cut -f1 -d"."); bedtools getfasta -fi ~/genome/Cobs.alpha.2.1.fa -bed $i -s |fold -w 60 > $base.fa;done

# add LG name to each sequence
for i in LG*.fa; do l=$(echo $i|cut -f1 -d"."); echo  $l $i; echo ">"$l > $l.final.fa; cat $i|grep -v ">" >> $l.final.fa;done

# merge them into one fasta file
cat LG*.final.fa|seqkit sort -l -r - > LGs.sorted.fa

# get sequences not placed in any LG
cat LG*bed |sort -V -k1,1 -k2,2 > regionsInLGs.bed

bedtools complement -i regionsInLGs.bed -g Cobs2.1.chrom-size.txt > regionsNotInLGs.bed

bedtools getfasta -fi ~/sciebo/My_PhD/Cobs_genomic_resources/genome/Cobs.alpha.2.1.fa -bed regionsNotInLGs.bed| seqkit sort -l -r - > regionsNotInLGs.sorted.fa

# merge LG groups and unplaced sequences
cat LGs.sorted.fa regionsNotInLGs.sorted.fa > Cobs3.0.tmp.fa


# And to remove sequences with only NNNs
seqtk comp Cobs3.0.tmp.fa|awk '{c=0;for(i=3;i<=6;++i){c+=$i}; gsub(">",""); if(c>0) print $1}'|xargs samtools faidx Cobs3.0.tmp.fa > Cobs3.0.clean.fa

picard NormalizeFasta I=Cobs3.0.clean.fa O=Cobs3.1.clean.fa

```
### Genome-wide recombination rates (gwRR)

Windowed recombination rates were estimated from the MSTmap output. We plotted genetic versus physical distance and visually checked correlations. Markers showing discrepencies between genetic and physical distances were removed. The final ordered markers file can be found in the datasets folder. I would like to thank Antonia Klein! Here code was a very great help here.

```R
# load libraries
library(tidyverse)
library(ggpubr)

# load data
lg <- read_delim("datasets/manuallyCurated.markers.forR.txt", col_names = F, delim = "\t")
colnames(lg) <- c("chr", "pos", "id", "distance1", "lg", "distance")
head (lg)

# reorder lg:
plotthis <- as.data.frame(lg[,c(5,6,3)])
colnames(plotthis) <- c("group", "position", "locus")
head(plotthis)

## by lg size
df <- plotthis %>% group_by(group) %>% summarise(max(position))
df <- df[order(df$`max(position)`, decreasing = T),]
df$group

## by no of markers in each lg
df2 <- plotthis %>% group_by(group) %>% summarise(length(position))
df2 <- df2[order(df2$`length(position)`, decreasing = T),]
df2$group

paste(shQuote(as.list(df$group)), collapse=", ")
lgs <- c('LG1', 'LG2', 'LG3', 'LG4', 'LG5', 'LG6', 'LG7', 'LG8', 'LG9', 'LG10', 'LG11',
         'LG12', 'LG13', 'LG14', 'LG15', 'LG16', 'LG17', 'LG18', 'LG19', 'LG20', 'LG21', 'LG22',
         'LG23', 'LG24', 'LG25', 'LG26', 'LG27', 'LG28', 'LG29')
paste(shQuote(as.list(sub("LG", "", lgs))), collapse=", ")

mcombsub <- reshape2::melt(lg, id=c("chr", "pos", "id", "lg", "distance1"))
mcombsub$CHR <- as.numeric(sub("LG", "", lg$lg))

mcombsub$CHR <- factor(mcombsub$CHR, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9',
                                                '10', '11', '12', '13', '14', '15', '16', '17',
                                                '18', '19', '20', '21', '22', '23', '24', '25', '26',
                                                '27', '28', '29'))

# plot physical vs genetic distances:
                                                                                ##########################################
                                                                                ### Plot physical vs genetic distances ###
                                                                                ##########################################


ggplot(filter(mcombsub, CHR %in% c('1', '2', '3', '4', '5', '6', '7', '8', '9',
                                   '10', '11', '12', '13', '14', '15', '16', '17',
                                   '18', '19', '20', '21', '22', '23', '24', '25', '26',
                                   '27', '28', '29')),
       aes(pos/1000000, value, color=variable))+  theme_bw()+
  geom_jitter(aes(), size=.4, alpha=.2)+
  geom_line(aes()) +
  stat_cor(method = "s",label.y.npc="top", label.x.npc = "left", cor.coef.name = "rho", cex=2, colour = "red")+
  scale_fill_manual(values=c("#66c2a5"))+
  scale_color_manual(values=c("#66c2a5"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
  )+
  scale_y_continuous(name=expression(paste("genetic distance (cM)", italic(" "))))+
  scale_x_continuous(name=expression(paste("physical distance (Mb)", italic(" "))))+
  theme(legend.position = "none")+
  # facet_grid(cols=vars(CHR),  scales = "free", space = "free_x")+
  facet_wrap(~CHR, ncol=5, nrow = 6, scales = "free")+
  theme(panel.spacing.x=unit(.5,"lines"),
        panel.spacing.y=unit(.5,"lines"),
        panel.background = element_blank(),
        axis.text.x = element_text(size=7, angle=45),
        #axis.ticks.x = element_blank(),
        panel.spacing = unit(.1, "cm"),
        #panel.background=element_rect(fill='white', colour='black'),
        strip.background=element_rect(fill='white', colour='black'))+ ggtitle("")+
  guides(colour = guide_legend(nrow = 1))

# dev.print(pdf,"~/sciebo/My_PhD/Cobs3.1/linkageMap/analyses/RRfinalMarkerSet/Figure.S2.pdf", height=10, width=8)

                                                                                ##################################################
                                                                                ### Plot Marker distribution along chromosomes ###
                                                                                ##################################################
# load data
lg <- read_delim("datasets/manuallyCurated.markers.forR.txt", col_names = F, delim = "\t")
colnames(lg) <- c("chr", "pos", "id", "distance1", "lg", "distance")

lg$scf = as.numeric(gsub("LG", "", lg$chr))

head (lg)

# load chromosome sizes
scf = read.table("datasets/Cobs3.1.chrom-size.txt")
scf = scf[seq(1,29,1),]
scf = cbind(scf, seq(1,29,1))
names(scf) = c("scf", "length", "number")

head(scf)
# plot

bp <- barplot(scf$length/1E6, border=NA, col="grey90",
              ylim=c(0,15), xlab="linkage group", ylab="length in Mb", names.arg=scf$number,
              cex.axis=.7, cex.names=.7, space = 10)

with(lg,
     segments(
       bp[scf,]-0.5,
       pos/1E6,
       bp[scf,]+0.5,
       pos/1E6,
       col="black",
       lwd=.5,
       lend="butt",
     )
)


                                                                                ##########################################
                                                                                ### Plot markers and linkage groups ###
                                                                                ##########################################
# load library
library("LinkageMapView")

# output file
outfile = file.path("sciebo/My_PhD/Cobs3.1/linkageMap/analyses/RRfinalMarkerSet/manuallyCurated.markers.forR.LGs.v2.pdf")

# plot
lmv.linkage.plot(plotthis,outfile, mapthese=c('LG7', 'LG2', 'LG3', 'LG5', 'LG1', 'LG8', 'LG21', 'LG9', 'LG15',
                                              'LG20', 'LG23', 'LG24', 'LG14', 'LG10', 'LG25', 'LG16', 'LG19',
                                              'LG11', 'LG4', 'LG27', 'LG13', 'LG22', 'LG18', 'LG29', 'LG12',
                                              'LG17', 'LG6', 'LG28', 'LG26'),
                 dupnbr=TRUE, lg.col = "#CCCCFF", lgperrow = 15, ruler = T, pdf.width = 56, pdf.height = 40)


                                                                                ######################################
                                                                                ### Computing gwRecombination Rates###
                                                                                ######################################

# recombination between recombining markers
rr1 <- lg %>%
  group_by(chr) %>%
  arrange(pos, .by_group = TRUE) %>%
  mutate(RR = distance - lag(distance, default = first(distance)),
         PhysicalDist = (pos - lag(pos, default = first(pos)))/1000000,
         pos2=lag(pos, default = first(pos)))

# cluster non recombining markers into regions
m1 <-
  lg %>%
  group_by(chr,distance) %>%
  arrange(pos, .by_group = TRUE) %>%
  filter(row_number() == 1 | row_number() == n())


# calculate recombination rate as: the ratio of the genetic to physical distance (in cM/Mb)
rr <- m1 %>%
  group_by(chr) %>%
  arrange(pos, .by_group = TRUE) %>%
  mutate(clusterGeneticDistMb = (distance - lag(distance, default = first(distance))),
         clusterPhysicalDistcM = (pos - lag(pos, default = first(pos)))/1000000,
         recombRatecMMb = (distance - lag(distance, default = first(distance)))/((pos - lag(pos, default = first(pos)))/1000000))

# save this!
write_delim(rr, "manuallyCurated.markers.RecomRateBetweenClustersGold.txt", delim = "\t")

# you need to edit the output in text editor to get recombination rates in a bed file and subsequently use bedtools intersect!
# bedtools makewindows -g chrom-size.txt -w 250000 -i srcwinnum |bedtools sort > 250kbwindows.bed
# bedtools intersect -a recombinationRateinBed.bed -b 250kbwindows.bed -wa -wb > recombinationRatein250kb.windows.txt


                                                                                ##################################################
                                                                                ### CALCULATE RECOMEBINATION IN SLIDING WINOWS ###
                                                                                ##################################################
data = read.table("recombinationRatein250kb.windows.txt", header=F)
head(data)

data = data[ ,c(1, 2, 3, 4, 6, 7, 8)]
names(data) <- c("scf","scf_position1","scf_position2", "RR",  "window_start", "window_end",  "window_name")

# calculate weighted mean of "RR" for each window:
##################################################
datalist = list()

for (j in unique(data$scf)) {
  print(j)
  data1 <- data[data$scf==j,]
  data1$window_name_short = as.numeric(gsub("^[^_]*_", "", data1$window_name))
  dt = data1
  dt$distance_bp <- data1$scf_position2-data1$scf_position1
  df = data.frame(window_name_short=seq(1,max(dt$window_name_short),1), window_size=NA, mean_RR=NA, weighted_mean_RR=NA)
  for(i in unique(dt$window_name_short)) {
    window = paste("window", i, sep="")
    assign(window, subset(dt, window_name_short==i))

    # if the window contains only one row, set window size to 1 and use its RR value for mean and weighted mean calculation
    if(nrow(get(paste("window",i,sep="")))==1) {
      df$window_size[i] <- "1"
      df$mean_RR[i] <- dt[dt$window_name_short==i, ]$RR
      df$weighted_mean_RR[i] <- dt[dt$window_name_short==i, ]$RR
    }

    else {
      df$window_size[i] <- nrow(get(paste("window",i,sep="")))

      # calculate mean RR value for the current window and store in 'df' data frame
      df$mean_RR[i] <- mean(dt[dt$window_name_short==i, ]$RR)

      # calculate weighted mean RR value
      # save rr of this window in a new vector:
      rr_window = paste("rr_window", i, sep="")
      assign(rr_window, dt[dt$window_name_short==i, ]$RR)

      # save physical distances of this window in a new vector:
      dist_window = paste("dist_window", i, sep="")
      assign(dist_window, dt[dt$window_name_short==i, ]$distance_bp)

      # calculate rr_window * dist_window
      product = paste("product", i, sep="")
      assign(product, get(paste("rr_window", i, sep="")) * get(paste("dist_window", i, sep="")))

      # add:
      sum = paste("sum", i, sep="")
      assign(sum, sum(get(paste("product", i, sep=""))))		

      # sum_distance:
      sum_distance = paste("sum_distance", i, sep="")
      assign(sum_distance, sum(get(paste("dist_window", i, sep=""))))		

      # weighted_mean:
      w_mean = paste("w_mean", i, sep="")
      assign(w_mean,  get(paste("sum", i, sep="")) /get(paste("sum_distance", i, sep="")) )

      df$weighted_mean_RR[i] <- get(paste("w_mean", i, sep=""))
    }
  }
  # merge data_seq:
  ##################

  merged = merge(data1, df, by="window_name_short")

  # add column with window mid position as midPos -values:
  merged$midPos = (as.numeric(merged$window_end) + as.numeric(merged$window_start))*1E-6/2

  merged_unique = merged[!duplicated(merged[,"window_name"]),]
  datalist[[j]] <- merged_unique

  # remove objects with names containing 'window', 'product', 'distance', 'sum', and 'mean'
  rm(list = ls(pattern = "window|product|distance|sum|mean"))  
}

windowedRR = do.call(rbind, datalist)

# add windows without markers:
nonoverlapping = read.table("250kbwindows.bed", header=F)
names(nonoverlapping) = c("scf","window_start","window_end","window_name")

df2 = merge(windowedRR, nonoverlapping, by="window_name", all.x=T, all.y=T)

# store this!
write_delim(df2, "finalRecombinationRateEstimates.250kbwindows.txt", delim = "\t")
```


### Population genomics using pool-seq data

For estimating nucleotide diversity (&#960;), Tajima's _D_ and genetic differentiation (_F_<sub>ST</sub>) from pool-seq data we followed workflow described before accessible [here](https://github.com/merrbii/CobsPopGenomics#iv--population-genomic-analyses). For nucleotide diversity at synonymous (&#960;<sub>S</sub>) and nonsynonymous (&#960;<sub>NS</sub>) sites, we used script, and available under [PoPoolation](https://sourceforge.net/projects/popoolation/) as follow. These metrics were calculated in nonoverlapping 250-kb windows.

```bash
# for Itabuna pool
popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 30 --gtf Cobs.alpha.v.3.1.geneannotation.1.5.polished.gtf --pileup A006200053_109203_S22_L002_CSFRD_indelsfree.pileup --output pool_Itabuna_syn-nonsyn_250kb250kb_nsl_p1.pi --window-size 250000 --step-size 250000 --fastq-type sanger --min-coverage 4 --max-coverage 75 --codon-table popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table popoolation_1.2.2/syn-nonsyn/nsl_p1.txt

# for Una pool
popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 50 --gtf Cobs.alpha.v.3.1.geneannotation.1.5.polished.gtf --pileup A006200053_109205_S23_L002_CSFRD_indelsfree.pileup --output pool_Una_syn-nonsyn_250kb250kb_nsl_p1.pi --window-size 250000 --step-size 250000 --fastq-type sanger --min-coverage 4 --max-coverage 77 --codon-table popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table popoolation_1.2.2/syn-nonsyn/nsl_p1.txt
```
### Molecular evolutionary rates

```bash
#####################################
# load genomes,CDSs and peptide sequences
# Solenopsis invicta
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/solenopsis_invicta/dna/Solenopsis_invicta.UNIL_Sinv_3.0.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/solenopsis_invicta/cds/Solenopsis_invicta.UNIL_Sinv_3.0.cds.all.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/solenopsis_invicta/pep/Solenopsis_invicta.UNIL_Sinv_3.0.pep.all.fa.gz

# monomorium pharaonis
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/monomorium_pharaonis_gca013373865v2/dna/Monomorium_pharaonis_gca013373865v2.GCA013373865v2.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/monomorium_pharaonis_gca013373865v2/cds/Monomorium_pharaonis_gca013373865v2.GCA013373865v2.cds.all.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/monomorium_pharaonis_gca013373865v2/pep/Monomorium_pharaonis_gca013373865v2.GCA013373865v2.pep.all.fa.gz

# ooceraea biroi
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/ooceraea_biroi_gca003672135v1/dna/Ooceraea_biroi_gca003672135v1.Obir_v5.4.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/ooceraea_biroi_gca003672135v1/cds/Ooceraea_biroi_gca003672135v1.Obir_v5.4.cds.all.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/ooceraea_biroi_gca003672135v1/pep/Ooceraea_biroi_gca003672135v1.Obir_v5.4.pep.all.fa.gz

# Nasonia vitripennis
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/nasonia_vitripennis/dna/Nasonia_vitripennis.Nvit_psr_1.1.dna.toplevel.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/nasonia_vitripennis/pep/Nasonia_vitripennis.Nvit_psr_1.1.pep.all.fa.gz
wget http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-54/metazoa/fasta/nasonia_vitripennis/cds/Nasonia_vitripennis.Nvit_psr_1.1.cds.all.fa.gz

# get the same for C. obscurior
ln -s ~/data/genomes/Cardiocondyla_obscurior/Cobs3.1/Cobs3.1.clean.fa .
ln -s ~/data/genomes/Cardiocondyla_obscurior/Cobs3.1/annotations/geneAnnotation/Cobs.alpha.v.3.1.geneannotation.1.5.polished.protein.fa .


#####################################
# use agat (https://github.com/NBISweden/AGAT) to get CDSs
agat_sp_extract_sequences.pl -g ~/data/genomes/Cardiocondyla_obscurior/Cobs3.1/annotations/geneAnnotation/Cobs.alpha.v.3.1.geneannotation.1.5.polished.gff3 -f ../genomes/Cobs3.1.clean.fa -t cds --output Cobs.alpha.v.3.1.geneannotation.1.5.polished.cds.fa

# Simplify headers by keeping only transcriptID:geneID, for instance >XM_011164587.3_gene:LOC105197943 or >XM_028189372.2_gene:LOC105840702 or >COBS16369-mRNA-1_gene:COBS16369
for i in Mpha Nvit Obir Sinv; do echo $i; awk -F' '  '{print (/^>/ ? $5"_"$4 : $0)}' $i.pep.fa|sed 's/transcript:/>/g' > tmp.txt && mv tmp.txt $i.pep.fa;done
cat pep/Cobs.pep.fa |awk -F' '  '{print (/^>/ ? $1"_"$2 : $0)}'|sed 's/=/:/g' > tmp.txt && mv tmp.txt pep/Cobs.pep.fa

# get longest isoform peptide sequence; speciesList.txt has species names (also used as file names), i.e. Mpha Nvit Obir Sinv Cobs
https://bioinf.shenwei.me/seqkit/
for i in $(cat speciesList.txt); do seqkit fx2tab -l pep/$i.pep.fa|sed 's/_gene:/\t/g' |sort -k2,2 -k4,4nr|sort -k2,2 -u -s|awk '{print ">"$1";gene="$2";length="$4"\n"$3}' > longestIso/longest.$i.pep.fa;done

# get longest isoform/transcript ID
for i in $(cat speciesList.txt );do grep "^>" longestIso/longest.$i.pep.fa|cut -f1 -d";"|sed 's/>//g' > cds/$i.longest_cds.txt ; done
for i in *.fa; do awk '/^>/ {$0=$1} 1' $i > tmp.fa && mv tmp.fa $i;done

# get longest isoform CDSs sequence
for i in $(cat ../speciesList.txt );do seqkit grep -n -f $i.longest_cds.txt $i.cds.fa > longest.$i.cds.fa ; done #removed -r to match full id

# unify headers
for i in $(cat speciesList.txt); do cat cds/longest.$i.cds.fa|sed "/^>/ s/$/_$i/" > cds/longest.simpleH.$i.cds.fa;done
for i in $(cat speciesList.txt); do cat pep/longest.$i.pep.fa|sed "s/;.*$/_$i/g" > pep/longest.simpleH.$i.pep.fa;done


#####################################
# prepare for OrthoFinder
mkdir orthofinder && mkdir inputs
ln -s ~/data/Cobs3.1/linkageMap/dNdS/inputs/pep/longest.simpleH.* .

# run OrthoFinder
nice orthofinder -f inputs/ -t 25 > orthofinder.log


#####################################
# prepare for sequence alignments
mkdir alignments

# get ortholgs prot
cat orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroups/Orthogroups_SingleCopyOrthologues.txt |while read group ; do cp orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroup_Sequences/"$group".fa alignments/"$group".faa ; done

# get corresponding cds
cat orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroups/Orthogroups_SingleCopyOrthologues.txt |while read group ; do cat speciesList.txt | while read species ; do grep "$species" ./alignments/"$group".faa | sed 's/>//g' | seqkit grep -n -f - cds/longest.simpleH."$species".cds.fa >> ./alignments/"$group".fna ; done ; done

# produce alignment
cat ~/data/Cobs3.1/linkageMap/dNdS/inputs/orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroups/Orthogroups_SingleCopyOrthologues.txt |nice parallel --jobs 25 --eta 'clustalo -i {1}.faa -o {1}.aln.faa'

# converts a alignments and corresponding cds sequences into a codon alignment
cat ~/data/Cobs3.1/linkageMap/dNdS/inputs/orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroups/Orthogroups_SingleCopyOrthologues.txt |time parallel --eta  pal2nal.pl {1}.aln.faa {1}.fna -output fasta -nogap ">" {1}.pal2nal.fasta "2>" {1}.pal2nal.fasta.log

# remove failed and erroneous alignments
grep "ERROR" *.log | sed 's/.pal2nal.*//g' > failed.txt
cat ~/data/Cobs3.1/linkageMap/dNdS/inputs/orthofinder/inputs/OrthoFinder/Results_Oct08_1/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | grep -v -f failed.txt > passedAlignments.txt

# filter codon alignments
ls *.fasta| nice parallel --jobs 7 --eta '/global/projects/programs/source/Gblocks_0.91b-new/Gblocks  {} -t=c -p=y' >> gblocks.stdout 2>&1

# run RAxML for each gene
ls *.fasta-gb|nice parallel --eta --jobs 30 'raxmlHPC -m GTRGAMMA -p 12345 -# 20 -s {} -n {.}.tree1 -T 1'
ls *.fasta-gb|nice parallel --eta --jobs 30 'raxmlHPC -m GTRGAMMA -p 12345 -b 12345 -# 100 -s {} -n {.}.tree2 -T 1'
ls *.fasta-gb|nice parallel --eta --jobs 30 'raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.{.}.tree1 -z RAxML_bootstrap.{.}.tree2 -n {.}.tree3 -T 1'

# Run aBSREL
cat passedAlignments.txt | while read group; do hyphy absrel --alignment sequence/"$group".pal2nal.fasta --tree bestTrees/RAxML_bipartitionsBranchLabels."$group".pal2nal.tree3 --output absrel/absrel-hyphy-"$group".json >> absrel/"$group".absrel.log 2>&1;done

# Run the FITMG94 model:
cat /home/m/merrbii/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/passedAlignments.txt |while read group; do nice ./HYPHYMP hyphy-analyses/FitMG94/FitMG94.bf --alignment ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/sequence/"$group".pal2nal.fasta --tree ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/bestTrees/RAxML_bipartitionsBranchLabels."$group".pal2nal.tree3 --output ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/FitMG94/fitMG94-hyphy-"$group".json --type local >> ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/FitMG94/"$group".fitMG94.log 2>&1;done
```

### Gene expression bias analysis

Here we used measures  previously calculated for RNA-seq data generated from worker, queen, wingless male, and winged male third instar larvae (n = 7 individuals each) ([Schrader et al., 2017](https://doi.org/10.1093/molbev/msw240)). But first we had to identify how genes in the current assembly/annotation (Cobs3.1) were referred to in the old assembly (Cobs1.4), used by Schrader et al.

```bash
# blast genes from Cobs1.4 against Cobs3.1
nice blastn -query Cobs3.1/Cobs3.1.genes.sequences.fa -subject Cobs1.4/Cobs1.4.genes.sequences.fa -outfmt "6 qseqid sseqid pident qlen length slen mismatch gapopen evalue bitscore" -out genesSeq.blast.results.txt

# get best uniq matches
cat genesSeq.blast.results.txt|awk -vOFS="\t" '{print $0,$5/$4}'|sort -k1,1 -k3,3n -k9,9n -k10,10nr -k11,11nr|sort -k1,1 -u -s > uniq.genesSeq.blast.results.txt

# get only those with at least 80 cov
awk '{if ($11 > 0.8) print $0}' uniq.genesSeq.blast.results.txt |cut -f1,2|tr -s ":"|tr ":" "\t"|cut -f1,4 > Cobs3.1.geneID.Cobs1.4.geneID.txt
```

The file obtained above contains two important columns; geneID(Cobs3.1) and geneid(Cobs1.4) + coordinates. Now in R merge this and the supplementary file generated by Schrader et al. and published along their paper.

```R
schrader <- read_delim("RNAseq/msw240_Supp/SupplementaryTable.csv", col_names = T, delim = ",")

mygenes <- read_delim("RNAseq/Cobs3.1.geneID.Cobs1.4.geneID.txt", col_names = F, delim = "\t")

colnames(mygenes) <- c("cobs31","chr", "start", "end","geneID")

merged <- merge(mygenes, schrader, by = "geneID", all.x = T)

## Add recombination rates:
# load rr
rr <- read_delim("windowedRR/finalRecombinationRateEstimates.250kbwindows.txt", col_names = T)
rr <- rr[,c(1,11)]
colnames(rr) <- c("window_name", "rr")


# load gene window names: a bed file containing 5 columns: 	 <geneID(Cobs3.1)> <chr> <start> <end> <window_name>
genewid <- read_delim("RNAseq/genesCobs3.1.matching.ExpressedgenesCobs1.4.RR.tmp.txt", col_names = F)
colnames(genewid) <- c("cobs31", "chr", "start", "end", "window_name")

# combine datasets
merged2 <- merge(genewid, rr, by = "window_name")
merged3 <- merge(merged, merged2, by = c("cobs31", "chr", "start", "end"))

# remove entries with no expression data:
merged3 <- merged3[!is.na(silver$`Expr(Ave)`),]

# remove duplicate genes: i.e. genes split in the new assembly
n_occur1 <- data.frame(table(merged3$geneID))
merged3 <- merged3[merged3$geneID %in% n_occur1$Var1[n_occur1$Freq < 2],]

data <- merged3

# store table for correlation analysis:
write.table(data,"RNAseq/finals/genesCobs3.1.matching.ExpressedgenesCobs1.4.Schraderetal.withRR.txt", quote = F, sep = "\t", row.names = F, col.names = T)

```
