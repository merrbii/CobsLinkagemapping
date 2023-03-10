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


* **_Get biallelic SNPs_**
```bash

vcftools --gzvcf gCobs3.1.linkageMapping.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out Cobs3.1.linkageMapping.snps

# compress and index
bgzip Cobs3.1.linkageMapping.snps.recode.vcf
tabix Cobs3.1.linkageMapping.snps.recode.vcf.gz
```
* **_Filter based on parental genotypes_**

Get heterozygous markers in F1 queen (F1Queen) that are homozygous fixed for different alleles in the parents or heterozygous in grandmother (PQueen)
```bash
gatk SelectVariants -V Cobs3.1.linkageMapping.snps.recode.vcf.gz -O biall.snps.noMendelianErrParents.vcf.gz -select 'vc.getGenotype("F1Queen").isHet() && (vc.getGenotype("PQueen").isHomRef() && vc.getGenotype("PMale").isHomVar() || vc.getGenotype("PQueen").isHet() && vc.getGenotype("PMale").isHom()) || vc.getGenotype("F1Queen").isHet() && (vc.getGenotype("PQueen").isHomVar() && vc.getGenotype("PMale").isHomRef() || vc.getGenotype("PQueen").isHet() && vc.getGenotype("PMale").isHom())' -select 'vc.getGenotype("F1Male").isHomVar() && (vc.getGenotype("PQueen").isHomVar() || vc.getGenotype("PQueen").isHet()) || vc.getGenotype("F1Male").isHomRef() && (vc.getGenotype("PQueen").isHomRef() || vc.getGenotype("PQueen").isHet())'
```
* **_Remove SNPs in problematic regions_**

In ants males are haploid and position called as heterozygous in these indicates a genotyping error and should be excluded. This is likely to happend in regions hard to align e.g. repeat rich regions. Get only positions where none of the males is heterozygous
```bash
gatk SelectVariants -V biall.snps.noMendelianErrParents.vcf.gz -O biall.snps.noMendelianErrParents.noHetMale.vcf.gz --invertSelect true -select 'vc.getGenotype("EM1").isHet()' -select 'vc.getGenotype("EM112").isHet()'  -select 'vc.getGenotype("EM2").isHet()' -select 'vc.getGenotype("EM3").isHet()' -select 'vc.getGenotype("EM35").isHet()' -select 'vc.getGenotype("EM36").isHet()'  -select 'vc.getGenotype("EM37").isHet()' -select 'vc.getGenotype("EM38").isHet()' -select 'vc.getGenotype("EM39").isHet()' -select 'vc.getGenotype("EM4").isHet()' -select 'vc.getGenotype("EM5").isHet()' -select 'vc.getGenotype("EM78").isHet()' -select 'vc.getGenotype("F1Male").isHet()' -select 'vc.getGenotype("PMale").isHet()'
```
* **_Discard genotypes showing mendelian errors_**

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
* **_Filter based on depth, missing data and allele frequency_**

Use [VCFtools](https://vcftools.github.io/index.html) to filter based on coverage, missingness and MAF.
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
Now switch to R and use [generate_phasedVCF.R](Rscripts/generate_phasedVCF.R) to combine annotations in `vcf.annotations.tsv` and the phased genotypes `phasedMarkers.noMissF1male.tsv`. After doing so, return to bash and produce the final clean vcf file with phased genotypes.
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

Windowed recombination rates were estimated from the MSTmap output. We plotted genetic versus physical distance and visually checked correlations. Markers showing discrepancies between genetic and physical distances were removed. The final ordered markers file can be found in the datasets folder. Use the [gwRecombinationRate.R](Rscripts/gwRecombinationRate.R) script to then calculate gwRR estimates as well as to visualize e.g. genetic versus physical distance correlations.


### Population genomics using pool-seq data

For estimating nucleotide diversity (&#960;), Tajima's _D_ and genetic differentiation (_F_<sub>ST</sub>) from pool-seq data we followed workflow described before accessible [here](https://github.com/merrbii/CobsPopGenomics#iv--population-genomic-analyses). For nucleotide diversity at synonymous (&#960;<sub>S</sub>) and nonsynonymous (&#960;<sub>NS</sub>) sites, we used script, and available under [PoPoolation](https://sourceforge.net/projects/popoolation/) as follow. These metrics were calculated in nonoverlapping 250-kb windows.

```bash
# for Itabuna pool
popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 30 --gtf Cobs.alpha.v.3.1.geneannotation.1.5.polished.gtf --pileup A006200053_109203_S22_L002_CSFRD_indelsfree.pileup --output pool_Itabuna_syn-nonsyn_250kb250kb_nsl_p1.pi --window-size 250000 --step-size 250000 --fastq-type sanger --min-coverage 4 --max-coverage 75 --codon-table popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table popoolation_1.2.2/syn-nonsyn/nsl_p1.txt

# for Una pool
popoolation_1.2.2/syn-nonsyn/Syn-nonsyn-sliding.pl --measure pi --pool-size 50 --gtf Cobs.alpha.v.3.1.geneannotation.1.5.polished.gtf --pileup A006200053_109205_S23_L002_CSFRD_indelsfree.pileup --output pool_Una_syn-nonsyn_250kb250kb_nsl_p1.pi --window-size 250000 --step-size 250000 --fastq-type sanger --min-coverage 4 --max-coverage 77 --codon-table popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table popoolation_1.2.2/syn-nonsyn/nsl_p1.txt
```
### Absolute divergence (_D_<sub>XY</sub>) calculation
For divergence analysis we used samples previously sequenced from four popolations; two representing the Old World lineage (collected from Taiwan and Leiden) and two of the New World lineage (Itabuna and Una) (see [Errbii et al. 2021](https://doi.org/10.1111/mec.16099)). Raw reads filtering, mapping and variant calling were performed as described previously ([here](https://github.com/merrbii/CobsPopGenomics)). We used GATK's GenotypeGVCFs with the `--all-sites` option to output invariant sites as well. _D_<sub>XY</sub> was then calculated using scripts written by **Dr. Simon H. Martin** that can be accessed [here](https://github.com/simonhmartin/genomics_general).


### Structural variants and short InDels analysis
For SVs and InDels we used sequencing data from 16 individual workers from five population of Cardiocondyla_obscurior_ previously generated (see [Errbii et al. 2021](https://doi.org/10.1111/mec.16099)). For InDels, we followed pipeline described before (further details [here](https://github.com/merrbii/CobsPopGenomics)) for raw reads filtering, mapping and variant calling. We then used GATK's _SelectVariants_ function with `--select-type-to-include INDEL` to only keep InDels. For filtering we used GATK's hard-filtering recommendations, `i.e. FS > 200.0, QD < 2.0, ReadPosRankSum < -20.0 and QUAL < 30.0`, and futher used VCFtools to keep bi-allelic ones `--max-alleles 2 and --min-alleles 2` with a MAF > 0.05 `--maf 0.05`.
For SVs, we followed steps detailed by **Dr. Claire M??rot** [here](https://github.com/clairemerot/SR_SV).

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
cat passedAlignments.txt | while read group; do hyphy absrel --alignment sequence/"$group".pal2nal.fasta-gb --tree bestTrees/RAxML_bipartitionsBranchLabels."$group".pal2nal.tree3 --output absrel/absrel-hyphy-"$group".json >> absrel/"$group".absrel.log 2>&1;done

# Run the FITMG94 model:
cat /home/m/merrbii/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/passedAlignments.txt |while read group; do nice ./HYPHYMP hyphy-analyses/FitMG94/FitMG94.bf --alignment ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/sequence/"$group".pal2nal.fasta-gb --tree ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/bestTrees/RAxML_bipartitionsBranchLabels."$group".pal2nal.tree3 --output ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/FitMG94/fitMG94-hyphy-"$group".json --type local >> ~/data/Cobs3.1/linkageMap/dNdS/inputs/hyphy/FitMG94/"$group".fitMG94.log 2>&1;done
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

The file obtained above contains two important columns; geneID(Cobs3.1) and geneid(Cobs1.4) + coordinates. Now in R merge this and the supplementary file generated by Schrader et al. and published along their paper. Use the [mergeRRexpressionbias.R](Rscripts/mergeRRexpressionbias.R) script to do so.
