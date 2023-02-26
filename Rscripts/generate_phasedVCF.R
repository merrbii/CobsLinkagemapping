# This short script merges annotations (ann) with phased genotypes (phased)
# generating a file that can be used to produce the final phased VCF


# load library
library(tidyverse)


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

