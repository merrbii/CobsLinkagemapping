# This script can be used to estimate windowed recombination rate along the genome based on MSTmap's output.
# I have attached the manually curated (i.e. after excluding ~200 markers showing discrepancies between physical/genetic distances)
# The file can be downloaded from the datasets folder.
# The first few sections of this script will produce useful plots, of which some are included as supplements to the manuscript
# I would like to thank Dr. Antonia Klein! Her code was a very great help here.

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
