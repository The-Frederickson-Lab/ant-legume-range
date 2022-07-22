Generalized mutualism promotes range expansion in both ant and plant
partners
================
Pooja Nathan and Megan Frederickson
Sys.Date()

# How does mutualism affect ant and plant range sizes?

This R Markdown document describes the dataset and code for:

Nathan P, Economo EP, Guenard B, Simonsen A, Frederickson ME.
Generalized mutualisms promote range expansion in both plant and ant
partners. In prep.

First, let’s load the packages we’ll use.

``` r
library(car)
library(ape)
library(geiger)
library(lme4)
library(nlme)
library(lmtest)
library(phytools)
library(caper)
library(tidyverse)
library(cowplot)
library(knitr)
library(taxize)
```

## Legume dataset

We obtained introduced and native range size data for legumes from
Simonsen et al. (2017).

Mutualistic trait data came from Weber et al. (2015) for extrafloral
nectaries (EFNs), Chomicki & Renner (2015) for domatia, Simonsen et
al. (2017) for nodulation, and Soudzilovskaia et al. (2020) for
mycorrhizae.

The next several chunks of code are slow to run, so they are not run
here, but they are included for reproducibility.

``` r
#Not run
#Legume range and nodulation data
range <- read.csv("legume_invasion_data_simonsenetal.csv") #Read in legume range data
names(range)[names(range) == "Species"] <- "Phy"
range$Phy <- as.character(gsub("_", " ", range$Phy))

#EFN data
EFN <- read.csv("EFNs_Weberatal_analysis_onlypresence.csv") #Read in EFN data
EFN$Phy <- as.character(gsub("_", " ", EFN$Phy))

#Domatia
domatia <- read.csv("domatia_chomickirenner_analysis_onlypresence.csv")  #Read in domatia data
names(domatia)[names(domatia) == 'ï..Phy'] <- 'Phy'
domatia$Phy <- as.character(gsub("_", " ", domatia$Phy))

#Mycorrhizae
mycorrhizae <- read.csv("mycorrhizae_sou.csv") #Read in legume species in Soudzilovskaia et al.'s mycorrhizae dataset
mycorrhizae$Phy <- as.character(gsub("_", " ", mycorrhizae$species))
mycorrhizae$In.range.dataset <- mycorrhizae$Phy %in% range$Phy #Check which taxa are in range dataset
mycorrhizae <- subset(mycorrhizae, In.range.dataset) #Subset to just taxa with mycorrhizal trait data that are also in legume range dataset
```

### Check taxonomy

We resolved misspelled names using the gnr_resolve function and checked
for synonyms using the synonyms function in the taxize package.

First, let’s check the taxonomy of the names of domatia-bearing plants.

Next, let’s check the taxonomy of EFN-bearing plants.

``` r
#Not run
#Do the same as above for the EFN data
EFN_resolve <- as.data.frame(gnr_resolve(EFN[, "Phy"], best_match_only=TRUE))
EFN_resolve$num_words <- str_count(EFN_resolve$matched_name, " ")+1
EFN <- merge(EFN, EFN_resolve, by.y = "user_supplied_name", by.x = "Phy", all.x=TRUE)
EFN$matched_name <- ifelse(EFN$num_words == 1, NA, EFN$matched_name)
EFN$matched_name <- ifelse(!is.na(EFN$matched_name), paste0(word(EFN$matched_name, 1, 1), " ", tolower(word(EFN$matched_name, 2, 2))), EFN$matched_name)
EFN$diff <- EFN$Phy == EFN$matched_name #Check changes between original and matched names

#Find synonyms for resolved names
#I had trouble getting this to work consistently so I retrieved the synonyms in batches
EFN_synonyms_pow1 <- synonyms(EFN[1:100,"matched_name"], db="pow") #Get synonyms
EFN_synonyms_pow2 <- synonyms(EFN[101:200,"matched_name"], db="pow")
EFN_synonyms_pow3 <- synonyms(EFN[201:400,"matched_name"], db="pow")
EFN_synonyms_pow4 <- synonyms(EFN[401:826,"matched_name"], db="pow")
EFN_synonyms_pow <- c(EFN_synonyms_pow1, EFN_synonyms_pow2, EFN_synonyms_pow3, EFN_synonyms_pow4)
EFN_syn <- do.call(rbind, EFN_synonyms_pow)
EFN_syn$matched_name <- gsub('[^-[:^punct:]]', '', (gsub('[[:digit:]]+', '', row.names(EFN_syn))), perl=TRUE)
EFN_syn <- subset(EFN_syn, rank == "SPECIES") #Keep on species (not varieties)
colnames(EFN_syn)[2] <- "synonym" #Fix column name

#Determine if original names and synonyms are in domatia, mycorrhizae, and range datasets
EFN_syn$syndomY <- EFN_syn$synonym %in% domatia$Phy 
EFN_syn$synrangeY <- EFN_syn$synonym %in% range$Phy
EFN_syn$synmycoY <- EFN_syn$synonym %in% mycorrhizae$Phy
write.csv(EFN_syn, "efn_synonyms.csv")
EFN_syn <- subset(EFN_syn, EFN_syn$syndomY | EFN_syn$synrangeY | EFN_syn$synmycoY) #Subset if synonyms match other datasets
EFN_syn$EFN <- 1 #Add trait
EFN_syn <- subset(EFN_syn, !is.na(synonym))
colnames(EFN_syn)[[2]] <- "Phy" #Fix column name

#Add synonyms to EFN dataset
EFN <- rbind(EFN[, c("EFN", "Phy", "matched_name")], EFN_syn[, c("EFN", "Phy","matched_name")]) #Merge
EFN$Phy <- trimws(EFN$Phy) #Trim white space from taxonomic names

#Check if each EFN name is in the EFN, mycorrhizae, and range datasets
EFN$RangeY <- EFN$Phy %in% range$Phy
EFN$domY <- EFN$Phy %in% domatia$Phy
EFN$mycoY <- EFN$Phy %in% mycorrhizae$Phy #No EFN-bearing plants in mycorrhizal dataset

write.csv(EFN, file="efn_resolved.csv", row.names=FALSE)
```

Merge the EFN and domatia datasets.

``` r
#Read in taxonomically resolved EFN and domatia datasets
EFN <- read.csv("efn_resolved.csv")
domatia <- read.csv("domatia_resolved.csv")

#Further cleaning of EFN and domatia datasets to 
colnames(EFN)[[3]] <- "matched_name_EFN"
colnames(domatia)[[3]] <- "matched_name_domatia"
domatia$diff <- domatia$Phy == domatia$matched_name_domatia 
EFN$diff <- EFN$Phy == EFN$matched_name_EFN
EFN <- subset(EFN, !is.na(diff))
EFN$matches <- EFN$RangeY+EFN$domY+EFN$mycoY
EFN <- subset(EFN, matches > 0)
EFN <- EFN[!duplicated(EFN), ]
EFN <- EFN[!duplicated(EFN$matched_name_EFN),]

#Merge 
EFN_dom <- merge(EFN[,c("EFN", "Phy", "matched_name_EFN")], domatia[,c("Domatia", "Phy")], by="Phy", all=TRUE)

#Add zeros instead of NAs for traits
EFN_dom$EFN <- ifelse(is.na(EFN_dom$EFN), 0, EFN_dom$EFN)
EFN_dom$Domatia <- ifelse(is.na(EFN_dom$Domatia), 0, EFN_dom$Domatia)
```

Next, we’ll check the taxonomy of the legumes for which we have data on
mycorrhizae.

``` r
#Not run
#Resolve names
mycorrhizae_resolve <- as.data.frame(gnr_resolve(mycorrhizae[, "Phy"], best_match_only=TRUE))
mycorrhizae_resolve$num_words <- str_count(mycorrhizae_resolve$matched_name, " ")+1
mycorrhizae <- merge(mycorrhizae, mycorrhizae_resolve, by.y = "user_supplied_name", by.x = "Phy", all.x=TRUE)
mycorrhizae$matched_name <- ifelse(mycorrhizae$num_words == 1, NA, mycorrhizae$matched_name)
mycorrhizae$matched_name <- ifelse(!is.na(mycorrhizae$matched_name), paste0(word(mycorrhizae$matched_name, 1, 1), " ", tolower(word(mycorrhizae$matched_name, 2, 2))), mycorrhizae$matched_name)

#Find synonyms for resolved names
mycorrhizae_synonyms_pow1 <- synonyms(unique(mycorrhizae$matched_name)[1:200], db="pow")
mycorrhizae_synonyms_pow2 <- synonyms(unique(mycorrhizae$matched_name)[201:400], db="pow")
mycorrhizae_synonyms_pow3 <- synonyms(unique(mycorrhizae$matched_name)[401:600], db="pow")
mycorrhizae_synonyms_pow4 <- synonyms(unique(mycorrhizae$matched_name)[601:800], db="pow")
mycorrhizae_syn <- append(mycorrhizae_synonyms_pow1, mycorrhizae_synonyms_pow2)
mycorrhizae_syn <- append(mycorrhizae_syn, mycorrhizae_synonyms_pow3)
mycorrhizae_syn <- append(mycorrhizae_syn, mycorrhizae_synonyms_pow4)

mycorrhizae_syn <- do.call(rbind, mycorrhizae_syn)
mycorrhizae_syn$matched_name <- gsub('[^-[:^punct:]]', '', (gsub('[[:digit:]]+', '', row.names(mycorrhizae_syn))), perl=TRUE)
mycorrhizae_syn <- subset(mycorrhizae_syn, rank == "SPECIES") #Keep on species (not varieties)
colnames(mycorrhizae_syn)[2] <- "synonym" #Fix column name

#Determine if original names and synonyms are in EFN, domatia, and range datasets
mycorrhizae_syn$syndomY <- mycorrhizae_syn$synonym %in% domatia$Phy 
mycorrhizae_syn$synrangeY <- mycorrhizae_syn$synonym %in% range$Phy
mycorrhizae_syn$synefnY <- mycorrhizae_syn$synonym %in% EFN$Phy
write.csv(mycorrhizae_syn, "mycorrhizae_synonyms.csv")

mycorrhizae_syn <- subset(mycorrhizae_syn, mycorrhizae_syn$syndomY | mycorrhizae_syn$synrangeY | mycorrhizae_syn$synefnY) #Subset if synonyms match other datasets
mycorrhizae_syn <- subset(mycorrhizae_syn, !is.na(synonym))
colnames(mycorrhizae_syn)[[2]] <- "Phy" #Fix column name

#Add synonyms to mycorrhizae dataset
mycorrhizae <- mycorrhizae[ ,-c(73:234)]
mycorrhizae_syn <- merge(mycorrhizae[, c("Phy", "mycorrhiza.type")], mycorrhizae_syn, by.x="Phy", by.y="matched_name", all.x=FALSE, all.y=TRUE)
colnames(mycorrhizae_syn)[[4]] <- "species"
mycorrhizae <- rbind(mycorrhizae[, c("mycorrhiza.type", "Phy", "species")], mycorrhizae_syn[, c("mycorrhiza.type", "Phy","species")]) #Merge
#mycorrhizae$Phy <- trimws(mycorrhizae$Phy) #Trim white space from taxonomic names

#Check if each mycorrhizae name is in the EFN, domatia, and range datasets
mycorrhizae$RangeY <- mycorrhizae$Phy %in% range$Phy
mycorrhizae$domY <- mycorrhizae$Phy %in% domatia$Phy
mycorrhizae$efnY <- mycorrhizae$Phy %in% EFN$Phy 

mycorrhizae$AM <- ifelse(mycorrhizae$mycorrhiza.type == "AM; others not addressed" | mycorrhizae$mycorrhiza.type == "AM; no others" | mycorrhizae$mycorrhiza.type == "EcM,AM", 1, 0)
mycorrhizae$EM <- ifelse(mycorrhizae$mycorrhiza.type == "EcM; others not addressed" | mycorrhizae$mycorrhiza.type == "EcM; no others" | mycorrhizae$mycorrhiza.type == "EcM,AM", 1, 0)

mycorrhizae.sum <- mycorrhizae %>% group_by(Phy, species, RangeY, domY, efnY) %>% summarize(n=n(), sum.AM = sum(AM), sum.EM=sum(EM))
mycorrhizae.sum$AM <- ifelse(mycorrhizae.sum$sum.AM > 0, "Y", "N")
mycorrhizae.sum$EM <- ifelse(mycorrhizae.sum$sum.EM > 0, "Y", "N")
mycorrhizae.sum$match <- ifelse(mycorrhizae.sum$Phy == mycorrhizae.sum$species, TRUE, FALSE)

write.csv(mycorrhizae.sum, file="mycorrhizae_resolved.csv", row.names=FALSE)
```

Merge the EFN and domatia dataset with the mycorrhizae dataset.

``` r
#Read in taxonomically resolved mycorrhizal dataset
mycorrhizae <- read.csv("mycorrhizae_resolved.csv")

#A little further cleaning of mycorrhizal dataset, to remove taxa that appear twice
mycorrhizae <- mycorrhizae %>% group_by(Phy, RangeY, domY, efnY) %>% summarize(n.records=sum(n), sum.AM = sum(sum.AM), sum.EM=sum(sum.EM))
mycorrhizae$AM <- ifelse(mycorrhizae$sum.AM > 0, "Y", "N")
mycorrhizae$EM <- ifelse(mycorrhizae$sum.EM > 0, "Y", "N")

#Merge 
EFN_dom_myco <- merge(mycorrhizae, EFN_dom, by="Phy", all=TRUE)
```

Finally, we’ll check the taxonomy of the range dataset.

``` r
#Not run
#Resolve names
range_resolve <- as.data.frame(gnr_resolve(range[, "Phy"], best_match_only=TRUE))
range_resolve$num_words <- str_count(range_resolve$matched_name, " ")+1
range <- merge(range, range_resolve, by.y = "user_supplied_name", by.x = "Phy", all.x=TRUE)
range$matched_name <- ifelse(range$num_words == 1, NA, range$matched_name)
range$matched_name <- ifelse(!is.na(range$matched_name), paste0(word(range$matched_name, 1, 1), " ", tolower(word(range$matched_name, 2, 2))), range$matched_name)
range$diff <- range$matched_name == range$Phy #Only one difference and spelling doesn't affect trait data (because taxon is not in EFN/dom/myco dataset)

write.csv(range, file="range_resolved.csv", row.names=FALSE)

#We don't really need to check synonyms, because we've already checked the synonyms of the matching datasets, and our goals with considering synonyms is just to make sure we match trait and range data correctly
```

Merge the mutualistic trait data with the range dataset.

``` r
range <- read.csv("range_resolved.csv")

legume_range_df <- merge(range, EFN_dom_myco, all.x="TRUE", all.y="FALSE", by= "Phy") #Put all the data in one dataframe
legume_range_df$EFN <- ifelse(is.na(legume_range_df$EFN), 0, legume_range_df$EFN) #Add zeros for NAs in EFN trait
legume_range_df$Domatia <- ifelse(is.na(legume_range_df$Domatia), 0, legume_range_df$Domatia) #Add zeros for NAs in domatia trait
write.csv(legume_range_df, file="legume_range_traits.csv", row.names = FALSE)
```

### Summarize legume dataset

We are finally finished cleaning the dataset, and now simply have to
summarize, visualize, and model.

First, let’s summarize how many species we have in each category. How
many legumes with vs. without EFNs do we have range size data for, and
how many introduced ranges have they been introduced to, on average?

``` r
#Make factors factors
legume_range_df$EFN <- as.factor(legume_range_df$EFN)
legume_range_df$Domatia <- as.factor(legume_range_df$Domatia)
legume_range_df$fixer <- as.factor(legume_range_df$fixer)
legume_range_df$AM <- as.factor(legume_range_df$AM)
legume_range_df$EM <- as.factor(legume_range_df$EM)
legume_range_df$annual <- as.numeric(legume_range_df$annual)

##Collapse all mycorrhizal fungi types into a single yes/no category
legume_range_df$myco <- ifelse(legume_range_df$AM == "Y" | legume_range_df$EM == "Y", 1, ifelse(legume_range_df$AM == "N" & legume_range_df$EM == "N", 0, NA))
legume_range_df$myco <- as.factor(legume_range_df$myco)

df <- legume_range_df

summary.efn <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(EFN) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.efn)
```

| EFN |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|-----:|--------------------:|------------------:|------------------:|
| 0   | 3695 |            1.366170 |          5.332059 |         0.0877178 |
| 1   |  280 |            5.257143 |         10.864345 |         0.6492688 |

How many legumes with vs. without domatia do we have range size data
for?

``` r
summary.dom <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(Domatia) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.dom)
```

| Domatia |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:--------|-----:|--------------------:|------------------:|------------------:|
| 0       | 3952 |           1.6467611 |          5.991534 |         0.0953080 |
| 1       |   23 |           0.5217391 |          1.201119 |         0.2504507 |

How many legumes that do vs. do not form nodules do we have range size
data for?

``` r
summary.fix <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(fixer) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.fix)
```

| fixer |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:------|-----:|--------------------:|------------------:|------------------:|
| 0     |  396 |            2.482323 |          6.839364 |         0.3436910 |
| 1     | 3579 |            1.547080 |          5.865714 |         0.0980483 |

How many legumes do vs. do not associate with mycorrhizae do we have
range size data for?

``` r
summary.myco <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(myco) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.myco)
```

| myco |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:-----|----:|--------------------:|------------------:|------------------:|
| 0    |  33 |            3.636364 |          7.192957 |         1.2521332 |
| 1    | 690 |            5.998551 |         11.970876 |         0.4557235 |

``` r
summary.myco2 <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(AM,EM) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.myco2)
```

| AM  | EM  |   n | mean_num_introduced | sd_num_introduced | se_num_introduced | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:----|:----|----:|--------------------:|------------------:|------------------:|---------------------:|-------------------:|-------------------:|
| N   | N   |  33 |           3.6363636 |          7.192957 |         1.2521332 |         3.194037e+12 |       8.008639e+12 |       1.394125e+12 |
| N   | Y   |  38 |           0.8684211 |          2.988004 |         0.4847182 |         4.655278e+11 |       1.604713e+12 |       2.603188e+11 |
| Y   | N   | 572 |           6.6311189 |         12.567518 |         0.5254743 |         6.482227e+12 |       1.410887e+13 |       5.899216e+11 |
| Y   | Y   |  80 |           3.9125000 |          9.169561 |         1.0251881 |         3.160244e+12 |       8.285881e+12 |       9.263897e+11 |

### Make figures

#### Number of introduced ranges

Our measure of ecological success for legumes is the number of new
ranges they have been successfully introduced to. We subsetted to
include only legumes with at least one introduced range. Each “inset”
shows the proportion of introduced species in with each trait.

##### Dots and whiskers

``` r
pt_size <- 3
y_limits <- c(-0.5, 15)
er_width <- 0.1
y_text <- -0.25
y_inset_limits <- c(0,1)

df$introducedY <- ifelse(df$num_introduced > 0, 1, 0) #Create binary variable for whether or not legume is introduced

#EFN figure
summary.efn <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(EFN) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.efn)
```

| EFN |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| 0   | 715 |             7.06014 |          10.33612 |         0.3865489 |
| 1   | 135 |            10.90370 |          13.55468 |         1.1666015 |

``` r
p_EFN <- ggplot(data=summary.efn, aes(x=EFN, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "***")

summary.efn2 <- ungroup(df %>% group_by(EFN, introducedY) %>% summarize(n=n()))
summary.efn2.wide <- spread(summary.efn2, key = introducedY, value=n)
colnames(summary.efn2.wide) <- c("EFN","Not_introduced",  "Introduced")
summary.efn2.wide$total <- summary.efn2.wide$Not_introduced+summary.efn2.wide$Introduced
summary.efn2.wide$prop.introduced <- summary.efn2.wide$Introduced/(summary.efn2.wide$total)
prop.efn <- paste0(summary.efn2.wide$Introduced, "/", summary.efn2.wide$total)

inset_p_EFN <- ggplot(data=summary.efn2.wide, aes(x=EFN, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=EFN, y=0.05, label=prop.efn), color="white")

#Domatia figure
summary.dom <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(Domatia) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.dom)
```

| Domatia |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:--------|----:|--------------------:|------------------:|------------------:|
| 0       | 844 |            7.710901 |          11.01974 |         0.3793152 |
| 1       |   6 |            2.000000 |           1.67332 |         0.6831301 |

``` r
p_dom <- ggplot(data=summary.dom, aes(x=Domatia, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "**")

summary.dom2 <- ungroup(df %>% group_by(Domatia, introducedY) %>% summarize(n=n()))
summary.dom2.wide <- spread(summary.dom2, key = introducedY, value=n)
colnames(summary.dom2.wide) <- c("Domatia","Not_introduced",  "Introduced")
summary.dom2.wide$total <- summary.dom2.wide$Not_introduced+summary.dom2.wide$Introduced
summary.dom2.wide$prop.introduced <- summary.dom2.wide$Introduced/summary.dom2.wide$total
prop.dom <- paste0(summary.dom2.wide$Introduced, "/", summary.dom2.wide$total)

inset_p_dom <- ggplot(data=summary.dom2.wide, aes(x=Domatia, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.dom), color="white")

#Nodules figure
summary.fix <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(fixer) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.fix)
```

| fixer |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:------|----:|--------------------:|------------------:|------------------:|
| 0     | 103 |            9.543689 |          10.63455 |         1.0478534 |
| 1     | 747 |            7.412316 |          11.02222 |         0.4032819 |

``` r
p_fix <- ggplot(data=summary.fix, aes(x=fixer, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "ns")

summary.fix2 <- ungroup(df %>% group_by(fixer, introducedY) %>% summarize(n=n()))
summary.fix2.wide <- spread(summary.fix2, key = introducedY, value=n)
colnames(summary.fix2.wide) <- c("fixer","Not_introduced",  "Introduced")
summary.fix2.wide$total <- summary.fix2.wide$Not_introduced+summary.fix2.wide$Introduced
summary.fix2.wide$prop.introduced <- summary.fix2.wide$Introduced/summary.fix2.wide$total
prop.fix <- paste0(summary.fix2.wide$Introduced, "/", summary.fix2.wide$total)
  
inset_p_fix <- ggplot(data=summary.fix2.wide, aes(x=fixer, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Nodules")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=fixer, y=0.05, label=prop.fix), color="white")

#Mycorrhizae figure
summary.AM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(AM) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.AM)
```

| AM  |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| N   |  22 |            6.954546 |          8.126836 |         1.7326471 |
| Y   | 308 |           13.331169 |         14.933754 |         0.8509296 |

``` r
p_AM <- ggplot(data=summary.AM, aes(x=AM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns")

summary.EM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(EM) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.EM)
```

| EM  |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| N   | 295 |           13.264407 |          14.92901 |         0.8692007 |
| Y   |  35 |            9.885714 |          11.88863 |         2.0095451 |

``` r
p_EM <- ggplot(data=summary.EM, aes(x=EM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("EM")+geom_text(aes(x=EM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns")

#Or
summary.AMEM2 <- ungroup(subset(df, !is.na(myco)) %>% group_by(AM, EM, introducedY) %>% summarize(n=n()))
summary.AMEM2.wide <- spread(summary.AMEM2, key = introducedY, value=n)
colnames(summary.AMEM2.wide) <- c("AM", "EM","Not_introduced",  "Introduced")
summary.AMEM2.wide$total <- summary.AMEM2.wide$Not_introduced+summary.AMEM2.wide$Introduced
summary.AMEM2.wide$prop.introduced <- summary.AMEM2.wide$Introduced/summary.AMEM2.wide$total
prop.AMEM <- paste0(summary.AMEM2.wide$Introduced, "/", summary.AMEM2.wide$total)

inset_p_AMEM <- ggplot(data=summary.AMEM2.wide, aes(x=AM, y=prop.introduced, fill=EM))+geom_bar(stat="identity", position=position_dodge())+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("AM")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1, y = 0.55, label = "*")+scale_fill_grey(labels=c("No", "Yes"))+geom_text(aes(x=c(0.77, 1.21, 1.78, 2.21), y=0.05, label=prop.AMEM), color="white")

fig1top <- plot_grid(inset_p_EFN, inset_p_dom, inset_p_fix, inset_p_AMEM, labels=c("AUTO"), nrow=1, rel_widths = c(1,1,1,2))
fig1bottom <-plot_grid(p_EFN, p_dom, p_fix, p_AM, p_EM, nrow=1, ncol=5, labels=c("E", "F", "G", "H", "I"))
fig1 <- plot_grid(fig1top, fig1bottom, nrow=2)
fig1
```

![](README_files/figure-gfm/Dot%20and%20whisker%20plots-1.png)<!-- -->

``` r
save_plot("Figure1.pdf", fig1, base_height=8, base_width=8)
```

#### Total introduced area

We might prefer to plot total introduced area, instead of the number of
introduced ranges.

``` r
pt_size <- 3
y_limits <- c(-500000, 7e+12)
er_width <- 0.1
y_text <- 0

summary.efn.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(EFN) %>% summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))

p_EFN_area <- ggplot(data=summary.efn.area, aes(x=EFN, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.dom.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(Domatia) %>% summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.dom.area)
```

| Domatia |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:--------|-----:|---------------------:|-------------------:|-------------------:|
| 0       | 3952 |         1.582757e+12 |       6.502918e+12 |       103442693450 |
| 1       |   23 |         4.312324e+11 |       1.017523e+12 |       212168306911 |

``` r
p_dom_area <- ggplot(data=summary.dom.area, aes(x=Domatia, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.fix.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(fixer) %>% summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.fix.area)
```

| fixer |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:------|-----:|---------------------:|-------------------:|-------------------:|
| 0     |  396 |         1.995301e+12 |       6.040271e+12 |       303535045685 |
| 1     | 3579 |         1.529711e+12 |       6.531601e+12 |       109178917397 |

``` r
p_fix_area <- ggplot(data=summary.fix.area, aes(x=fixer, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.myco.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(myco) %>% summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.myco.area)
```

| myco |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:-----|-----:|---------------------:|-------------------:|-------------------:|
| 0    |   33 |         3.194037e+12 |       8.008639e+12 |       1.394125e+12 |
| 1    |  690 |         5.765715e+12 |       1.325680e+13 |       5.046778e+11 |
| NA   | 3252 |         6.707344e+11 |       2.993384e+12 |       5.249128e+10 |

``` r
p_myco_area <- ggplot(data=subset(summary.myco.area, !is.na(myco)), aes(x=myco, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=myco, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("Mycorrhizae")+geom_text(aes(x=myco, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

fig1_area <- plot_grid(p_EFN_area, p_dom_area, p_myco_area, p_fix_area, nrow=1, labels="AUTO")
fig1_area
```

![](README_files/figure-gfm/Dot%20and%20whisker%20plots%20for%20introduced%20area-1.png)<!-- -->

``` r
save_plot("Figure1alt.pdf", fig1_area, base_height=4, base_width=8)
```

#### Total native area

We can also plot the effect of the same mutualistic traits on the native
range size of legumes.

``` r
pt_size <- 3
y_limits <- c(-500000, 13e+12)
er_width <- 0.1
y_text <- 0

summary.efn.dom.native <- ungroup(df %>% group_by(EFN, Domatia) %>% summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_EFN_dom_native <- ggplot(data=summary.efn.dom.native, aes(x=EFN, y=mean_native, color=Domatia))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native, color=Domatia), width=er_width)+ geom_line(aes(group=Domatia), linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("EFNs")+geom_text(aes(x=c(0.8,1.2,1.8,2.2), y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+scale_color_grey(labels=c("No", "Yes"))+theme(legend.position = c(0.1, 0.8))+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")+annotate("text", x = 1.5, y = 3.9e+12, label = "***")

summary.dom.native <- ungroup(df %>% group_by(Domatia) %>% summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_dom_native <- ggplot(data=summary.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.fix.native <- ungroup(df %>% group_by(fixer) %>% summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_fix_native <- ggplot(data=summary.fix.native, aes(x=fixer, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

summary.AM.native <- ungroup(df %>% group_by(AM) %>% summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.AM.native)
```

| AM  |    n |  mean_native |    sd_native |    se_native |
|:----|-----:|-------------:|-------------:|-------------:|
| N   |   71 | 6.934794e+12 | 7.655628e+12 | 908555938765 |
| Y   |  652 | 1.025923e+13 | 1.045558e+13 | 409472087225 |
| NA  | 3252 | 4.766370e+12 | 5.005758e+12 |  87779793963 |

``` r
p_AM_native <- ggplot(data=subset(summary.AM.native, !is.na(AM)), aes(x=AM, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

summary.EM.native <- ungroup(df %>% group_by(EM) %>% summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.EM.native)
```

| EM  |    n |  mean_native |    sd_native |    se_native |
|:----|-----:|-------------:|-------------:|-------------:|
| N   |  605 | 1.090616e+13 | 1.083941e+13 | 440684520332 |
| Y   |  118 | 5.058614e+12 | 3.946714e+12 | 363324484094 |
| NA  | 3252 | 4.766370e+12 | 5.005758e+12 |  87779793963 |

``` r
p_EM_native <- ggplot(data=subset(summary.EM.native, !is.na(EM)), aes(x=EM, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=EM, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("EM")+geom_text(aes(x=EM, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

fig2 <- plot_grid(p_EFN_dom_native, p_fix_native, p_AM_native, p_EM_native, nrow=1, labels="AUTO", rel_widths = c(1.8, 1, 1, 1))
fig2
```

![](README_files/figure-gfm/Dots%20and%20whiskers%20plots%20for%20native%20area-1.png)<!-- -->

``` r
save_plot("Figure2.pdf", fig2, base_height = 4, base_width = 12)
```

### Statistical models

#### Preparing dataset for pgls analysis

``` r
df$Phy2 <- paste0(as.character(word(legume_range_df$Phy, 1, 1)), "_", as.character(word(legume_range_df$Phy, 2, 2)))
zanne <- read.tree("Vascular_Plants_rooted.dated.tre") #reading in Zanne et al. 2014 plant phylogeny
phyint <- intersect(zanne$tip.label, df$Phy2)  
phydiff <- setdiff(zanne$tip.label, df$Phy2)
pruned.tree.pgls <- drop.tip(zanne, phydiff) #dropping tips not in the dataset

range_pgls <- df[df$Phy2 %in% phyint, ]
#colnames(range_pgls)
#which(colSums(is.na(range_pgls))>0) #Check which columns have NAs
range_pgls <-range_pgls[,-c(16:23, 25, 27)] #Remoove some unneeded columns
range_pgls <- range_pgls[complete.cases(range_pgls), ] #removing NA elements
```

#### PGLS models for EFN, domatia, and rhizobia

The problem with PGLS models is that we lose a lot of trait data for
species not in the phylogeny. We either need a better phylogeny or a
non-phylogenetic model.

``` r
#PGLS of number of introduced ranges as response variable with EFNs, domatia, and nodules, 
#interaction and other covariates
#Interactions that were not significant were removed

pgls1 <- gls(log(num_introduced + 1) ~  EFN+fixer+Domatia+total_area_native + abs_lat_native + uses_num_uses + annual, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method = "ML", data = range_pgls) 
summary(pgls1)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(num_introduced + 1) ~ EFN + fixer + Domatia + total_area_native +      abs_lat_native + uses_num_uses + annual 
    ##   Data: range_pgls 
    ##        AIC      BIC    logLik
    ##   2677.736 2728.794 -1328.868
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.2880877 
    ## 
    ## Coefficients:
    ##                        Value  Std.Error  t-value p-value
    ## (Intercept)        0.2328939 0.18310985  1.27188  0.2037
    ## EFN1               0.4666291 0.07397529  6.30790  0.0000
    ## fixer1             0.1029162 0.11425070  0.90079  0.3679
    ## Domatia1          -0.2729300 0.31038817 -0.87932  0.3794
    ## total_area_native  0.0000000 0.00000000 -1.89313  0.0586
    ## abs_lat_native    -0.0027053 0.00203418 -1.32990  0.1838
    ## uses_num_uses      0.4180435 0.01243464 33.61928  0.0000
    ## annual             0.0565769 0.06690117  0.84568  0.3979
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN1   fixer1 Domat1 ttl_r_ abs_l_ uss_n_
    ## EFN1               0.005                                          
    ## fixer1            -0.244 -0.069                                   
    ## Domatia1           0.002  0.000 -0.097                            
    ## total_area_native -0.065 -0.021 -0.041 -0.005                     
    ## abs_lat_native    -0.216  0.024 -0.058  0.054  0.081              
    ## uses_num_uses     -0.089 -0.119  0.071  0.020 -0.424  0.030       
    ## annual            -0.039 -0.008 -0.025  0.008  0.080  0.091  0.046
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.5299088 -0.3834094 -0.2998552  0.3561601  3.2716673 
    ## 
    ## Residual standard error: 0.810807 
    ## Degrees of freedom: 1219 total; 1211 residual

``` r
#Diagnostic plots
plot(pgls1$residuals, pgls1$fitted)
```

![](README_files/figure-gfm/Legume%20pgls%20models-1.png)<!-- -->

``` r
qqnorm(pgls1$residuals)
qqline(pgls1$residuals)
```

![](README_files/figure-gfm/Legume%20pgls%20models-2.png)<!-- -->

``` r
#PGLS of total introduced area
pgls2 <- gls(log(total_area_introduced/1e+6 + 1) ~  EFN+Domatia+fixer+
            total_area_native + abs_lat_native + uses_num_uses+annual, 
            correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), 
            method = "ML", data = range_pgls) 
summary(pgls2)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(total_area_introduced/1e+06 + 1) ~ EFN + Domatia + fixer +      total_area_native + abs_lat_native + uses_num_uses + annual 
    ##   Data: range_pgls 
    ##        AIC      BIC    logLik
    ##   7635.715 7686.773 -3807.858
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.1901358 
    ## 
    ## Coefficients:
    ##                       Value Std.Error   t-value p-value
    ## (Intercept)        1.765915 1.1694960  1.509979  0.1313
    ## EFN1               3.176319 0.5673187  5.598827  0.0000
    ## Domatia1          -1.288304 2.3447548 -0.549441  0.5828
    ## fixer1             0.656741 0.8164576  0.804378  0.4213
    ## total_area_native  0.000000 0.0000000  0.260467  0.7945
    ## abs_lat_native     0.021659 0.0151153  1.432890  0.1521
    ## uses_num_uses      2.252582 0.0949824 23.715791  0.0000
    ## annual             0.834530 0.5069873  1.646056  0.1000
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN1   Domat1 fixer1 ttl_r_ abs_l_ uss_n_
    ## EFN1               0.006                                          
    ## Domatia1           0.000  0.003                                   
    ## fixer1            -0.297 -0.071 -0.096                            
    ## total_area_native -0.077 -0.025 -0.004 -0.045                     
    ## abs_lat_native    -0.247  0.026  0.061 -0.072  0.085              
    ## uses_num_uses     -0.107 -0.117  0.021  0.081 -0.425  0.029       
    ## annual            -0.045 -0.009  0.010 -0.030  0.073  0.087  0.048
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.8666296 -0.6526755 -0.4946179  0.7096023  2.3323612 
    ## 
    ## Residual standard error: 5.901819 
    ## Degrees of freedom: 1219 total; 1211 residual

``` r
#Diagnostic plots
plot(pgls2$residuals, pgls2$fitted)
```

![](README_files/figure-gfm/Legume%20pgls%20models-3.png)<!-- -->

``` r
qqnorm(pgls2$residuals)
qqline(pgls2$residuals)
```

![](README_files/figure-gfm/Legume%20pgls%20models-4.png)<!-- -->

``` r
#Repeating for native area
#PGLS with both EFN and domatia presence and fixer, interaction and covariates
pgls3 <- gls(log((total_area_native/1e+6) + 1) ~ EFN*Domatia + fixer+ abs_lat_native+annual + uses_num_uses, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method = "ML", data = range_pgls)
summary(pgls3)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log((total_area_native/1e+06) + 1) ~ EFN * Domatia + fixer +      abs_lat_native + annual + uses_num_uses 
    ##   Data: range_pgls 
    ##       AIC      BIC   logLik
    ##   4099.36 4150.418 -2039.68
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.4868966 
    ## 
    ## Coefficients:
    ##                    Value Std.Error  t-value p-value
    ## (Intercept)    14.872340 0.4448421 33.43285  0.0000
    ## EFN1           -0.084274 0.1311168 -0.64274  0.5205
    ## Domatia1        0.520235 0.6360372  0.81793  0.4136
    ## fixer1          0.052313 0.2243069  0.23322  0.8156
    ## abs_lat_native -0.007530 0.0037735 -1.99560  0.0462
    ## annual         -0.042524 0.1202243 -0.35370  0.7236
    ## uses_num_uses   0.247038 0.0201623 12.25247  0.0000
    ## EFN1:Domatia1  -3.329667 1.3613747 -2.44581  0.0146
    ## 
    ##  Correlation: 
    ##                (Intr) EFN1   Domat1 fixer1 abs_l_ annual uss_n_
    ## EFN1            0.003                                          
    ## Domatia1        0.004  0.032                                   
    ## fixer1         -0.180 -0.069 -0.094                            
    ## abs_lat_native -0.165  0.026  0.040 -0.043                     
    ## annual         -0.026 -0.003  0.005 -0.016  0.087              
    ## uses_num_uses  -0.093 -0.146 -0.001  0.047  0.070  0.092       
    ## EFN1:Domatia1  -0.004 -0.080 -0.458  0.038 -0.004  0.001  0.030
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -7.3357948 -0.3211807  0.1350769  0.6054320  1.6110388 
    ## 
    ## Residual standard error: 1.646054 
    ## Degrees of freedom: 1219 total; 1211 residual

``` r
View(subset(range_pgls, EFN == 1 & Domatia == 1))

#Diagnostic plots
plot(pgls3$residuals, pgls3$fitted)
```

![](README_files/figure-gfm/Legume%20pgls%20models-5.png)<!-- -->

``` r
qqnorm(pgls3$residuals)
qqline(pgls3$residuals)
```

![](README_files/figure-gfm/Legume%20pgls%20models-6.png)<!-- -->

#### PGLS models for mycorrhizae

``` r
range_myco <- subset(df, !is.na(myco))
range_myco$Phy2 <- as.character(range_myco$Phy2)
phyint1 <- intersect(zanne$tip.label, range_myco$Phy2)  
phydiff1 <- setdiff(zanne$tip.label, range_myco$Phy2)
pruned.myco.pgls <- drop.tip(zanne, phydiff1) #dropping tips not in the dataset

range_myco_pgls <- range_myco[range_myco$Phy2 %in% phyint1, ]
which(colSums(is.na(range_myco_pgls))>0) #Check which columns have NAs
```

    ## total_area_native        lat_native    abs_lat_native  matched_name_EFN 
    ##                 5                 6                 7                25

``` r
range_myco_pgls <-range_myco_pgls[,-c(25)] #Remove some unneeded columns
range_myco_pgls <- range_myco_pgls[complete.cases(range_myco_pgls), ] #removing NA elements

#PGLS of number of introduced ranges as response variable mycorrhizae and covariates as predictors
pgls4 <- gls(log(num_introduced + 1) ~ AM+EM + total_area_native + abs_lat_native + annual + uses_num_uses, correlation = corPagel(1, phy = pruned.myco.pgls, form = ~ Phy2), method = "ML", data = range_myco_pgls) 
summary(pgls4)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(num_introduced + 1) ~ AM + EM + total_area_native + abs_lat_native +      annual + uses_num_uses 
    ##   Data: range_myco_pgls 
    ##        AIC      BIC    logLik
    ##   1020.099 1055.655 -501.0494
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.3627128 
    ## 
    ## Coefficients:
    ##                        Value  Std.Error   t-value p-value
    ## (Intercept)        0.5736411 0.31281273  1.833816  0.0675
    ## AMY                0.0262398 0.17302781  0.151651  0.8795
    ## EMY                0.1225168 0.15126770  0.809933  0.4185
    ## total_area_native  0.0000000 0.00000000 -2.520128  0.0121
    ## abs_lat_native    -0.0092471 0.00417507 -2.214838  0.0274
    ## annual             0.2869247 0.17870901  1.605541  0.1092
    ## uses_num_uses      0.4207177 0.02080559 20.221377  0.0000
    ## 
    ##  Correlation: 
    ##                   (Intr) AMY    EMY    ttl_r_ abs_l_ annual
    ## AMY               -0.478                                   
    ## EMY               -0.147  0.215                            
    ## total_area_native -0.104  0.001  0.060                     
    ## abs_lat_native    -0.280  0.027 -0.136 -0.002              
    ## annual            -0.077  0.033  0.019  0.199  0.119       
    ## uses_num_uses     -0.090 -0.104 -0.008 -0.295  0.148 -0.048
    ## 
    ## Standardized residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.69416567 -0.48205754 -0.08573794  0.66570930  2.67010580 
    ## 
    ## Residual standard error: 1.019246 
    ## Degrees of freedom: 384 total; 377 residual

#### Mixed models

Another possible approach is to use mixed models with legume genus as a
random effect to account for the non-independence of species in a genus.
The response variables (number of introduced ranges and total introduced
area) are very non-normal so we fit two models: a binomial GLMM
modelling whether or not a legume species is introduced, and linear
mixed models of the number of introduced ranges or total introduced area
for introduced species only.

``` r
df$genus <- as.factor(word(df$Phy2, 1, 1, sep="_")) #Extract legume genus

#Fit binomial model for whether or not a legume species has been introduced
binomial1 <- glmer(introducedY~EFN+Domatia+fixer+scale(total_area_native)+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=df, family="binomial")
summary(binomial1) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: introducedY ~ EFN + Domatia + fixer + scale(total_area_native) +  
    ##     annual + scale(abs_lat_native) + uses_num_uses + (1 | genus)
    ##    Data: df
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2529.6   2585.1  -1255.8   2511.6     3523 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -12.3796  -0.3510  -0.2260  -0.1257   6.3465 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 1.233    1.111   
    ## Number of obs: 3532, groups:  genus, 440
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -3.13006    0.26259 -11.920  < 2e-16 ***
    ## EFN1                      1.23977    0.20842   5.948 2.71e-09 ***
    ## Domatia1                  0.40484    0.81073   0.499  0.61753    
    ## fixer1                    0.45982    0.25121   1.830  0.06719 .  
    ## scale(total_area_native)  0.17357    0.06160   2.818  0.00483 ** 
    ## annual                    0.40501    0.17470   2.318  0.02043 *  
    ## scale(abs_lat_native)     0.07839    0.06955   1.127  0.25971    
    ## uses_num_uses             1.06957    0.05478  19.525  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 scl(t__) annual scl(b__)
    ## EFN1        -0.065                                              
    ## Domatia1    -0.023 -0.020                                       
    ## fixer1      -0.867  0.017 -0.006                                
    ## scl(ttl_r_)  0.137 -0.027  0.006 -0.086                         
    ## annual      -0.037 -0.002  0.013 -0.064  0.022                  
    ## scl(bs_lt_)  0.170  0.046  0.078 -0.173  0.092    0.034         
    ## uses_num_ss -0.361  0.011  0.032  0.129 -0.336    0.116  0.004

``` r
Anova(binomial1, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)              142.0879  1  < 2.2e-16 ***
    ## EFN                       35.3832  1  2.708e-09 ***
    ## Domatia                    0.2494  1   0.617527    
    ## fixer                      3.3504  1   0.067187 .  
    ## scale(total_area_native)   7.9404  1   0.004834 ** 
    ## annual                     5.3747  1   0.020431 *  
    ## scale(abs_lat_native)      1.2703  1   0.259706    
    ## uses_num_uses            381.2143  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial1)
```

![](README_files/figure-gfm/Legume%20mixed%20models-1.png)<!-- -->

``` r
#Fit binomial model without random effect
binomial2 <- glm(introducedY~EFN+Domatia+fixer+scale(total_area_native)+annual+scale(abs_lat_native)+uses_num_uses, data=df, family="binomial")
summary(binomial2) 
```

    ## 
    ## Call:
    ## glm(formula = introducedY ~ EFN + Domatia + fixer + scale(total_area_native) + 
    ##     annual + scale(abs_lat_native) + uses_num_uses, family = "binomial", 
    ##     data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.1562  -0.5233  -0.4089  -0.3784   2.3544  
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.47854    0.18642 -13.295  < 2e-16 ***
    ## EFN1                      1.38962    0.16999   8.175 2.96e-16 ***
    ## Domatia1                  0.21840    0.74260   0.294  0.76868    
    ## fixer1                    0.10345    0.18406   0.562  0.57409    
    ## scale(total_area_native)  0.17605    0.05358   3.286  0.00102 ** 
    ## annual                    0.63910    0.13610   4.696 2.66e-06 ***
    ## scale(abs_lat_native)     0.07651    0.04981   1.536  0.12455    
    ## uses_num_uses             0.99740    0.04689  21.270  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3811.0  on 3531  degrees of freedom
    ## Residual deviance: 2653.2  on 3524  degrees of freedom
    ##   (443 observations deleted due to missingness)
    ## AIC: 2669.2
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
Anova(binomial2, type=3)
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: introducedY
    ##                          LR Chisq Df Pr(>Chisq)    
    ## EFN                         62.89  1  2.190e-15 ***
    ## Domatia                      0.08  1  0.7728068    
    ## fixer                        0.32  1  0.5722022    
    ## scale(total_area_native)    10.92  1  0.0009493 ***
    ## annual                      20.94  1  4.743e-06 ***
    ## scale(abs_lat_native)        2.36  1  0.1243614    
    ## uses_num_uses              724.56  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-2.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-3.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-4.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-5.png)<!-- -->

``` r
#Compare models
anova(binomial1, binomial2) #Random effect improves model fit to data
```

    ## Data: df
    ## Models:
    ## binomial2: introducedY ~ EFN + Domatia + fixer + scale(total_area_native) + annual + scale(abs_lat_native) + uses_num_uses
    ## binomial1: introducedY ~ EFN + Domatia + fixer + scale(total_area_native) + annual + scale(abs_lat_native) + uses_num_uses + (1 | genus)
    ##           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
    ## binomial2    8 2669.2 2718.6 -1326.6   2653.2                         
    ## binomial1    9 2529.6 2585.2 -1255.8   2511.6 141.58  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Fit a linear mixed model for how many new ranges an introduced legume has established in, and how much total area they cover
legume_range_df_introducedY <- subset(df, num_introduced >0) #Filter to species with 1+ introduced ranges

#Number of introduced ranges
lmer1 <- lmer(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+uses_num_uses+(1|genus), data=legume_range_df_introducedY)
summary(lmer1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) +  
    ##     scale(total_area_native) + annual + uses_num_uses + (1 |      genus)
    ##    Data: legume_range_df_introducedY
    ## 
    ## REML criterion at convergence: 2200.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.88825 -0.73180 -0.03294  0.67398  2.47609 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.1172   0.3423  
    ##  Residual             0.7836   0.8852  
    ## Number of obs: 813, groups:  genus, 195
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.61854    0.13344   4.635
    ## EFN1                      0.38845    0.10134   3.833
    ## Domatia1                 -1.08782    0.40985  -2.654
    ## fixer1                   -0.13400    0.13099  -1.023
    ## scale(abs_lat_native)    -0.03533    0.04029  -0.877
    ## scale(total_area_native) -0.02947    0.03668  -0.803
    ## annual                    0.06725    0.10864   0.619
    ## uses_num_uses             0.32059    0.01644  19.502
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 scl(b__) scl(t__) annual
    ## EFN1        -0.086                                              
    ## Domatia1     0.028 -0.024                                       
    ## fixer1      -0.876  0.025 -0.025                                
    ## scl(bs_lt_)  0.105  0.109  0.063 -0.115                         
    ## scl(ttl_r_)  0.219 -0.022  0.005 -0.130  0.031                  
    ## annual      -0.038  0.043 -0.005 -0.080 -0.088    0.127         
    ## uses_num_ss -0.388 -0.064 -0.041  0.100  0.044   -0.379    0.050

``` r
Anova(lmer1, type=3,)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               21.4862  1  3.564e-06 ***
    ## EFN                       14.6936  1  0.0001265 ***
    ## Domatia                    7.0448  1  0.0079494 ** 
    ## fixer                      1.0465  1  0.3063054    
    ## scale(abs_lat_native)      0.7689  1  0.3805648    
    ## scale(total_area_native)   0.6454  1  0.4217503    
    ## annual                     0.3832  1  0.5358805    
    ## uses_num_uses            380.3192  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer1)
```

![](README_files/figure-gfm/Legume%20mixed%20models-6.png)<!-- -->

``` r
#Number of introduced ranges without random effect of genus
lm2 <- lm(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+uses_num_uses, data=legume_range_df_introducedY)
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + 
    ##     scale(total_area_native) + annual + uses_num_uses, data = legume_range_df_introducedY)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -2.98319 -0.73804 -0.06007  0.62360  2.50352 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               0.71351    0.11200   6.371 3.17e-10 ***
    ## EFN1                      0.33522    0.09148   3.664 0.000264 ***
    ## Domatia1                 -1.19195    0.42215  -2.824 0.004867 ** 
    ## fixer1                   -0.21908    0.10628  -2.061 0.039594 *  
    ## scale(abs_lat_native)    -0.02295    0.03505  -0.655 0.512783    
    ## scale(total_area_native) -0.02186    0.03558  -0.615 0.539010    
    ## annual                    0.26267    0.09714   2.704 0.006993 ** 
    ## uses_num_uses             0.31893    0.01651  19.321  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9347 on 805 degrees of freedom
    ##   (37 observations deleted due to missingness)
    ## Multiple R-squared:  0.3685, Adjusted R-squared:  0.363 
    ## F-statistic: 67.11 on 7 and 805 DF,  p-value: < 2.2e-16

``` r
Anova(lm2, type=3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(num_introduced)
    ##                          Sum Sq  Df  F value    Pr(>F)    
    ## (Intercept)               35.46   1  40.5834 3.167e-10 ***
    ## EFN                       11.73   1  13.4282 0.0002641 ***
    ## Domatia                    6.97   1   7.9722 0.0048673 ** 
    ## fixer                      3.71   1   4.2490 0.0395937 *  
    ## scale(abs_lat_native)      0.37   1   0.4288 0.5127835    
    ## scale(total_area_native)   0.33   1   0.3777 0.5390101    
    ## annual                     6.39   1   7.3123 0.0069930 ** 
    ## uses_num_uses            326.15   1 373.2889 < 2.2e-16 ***
    ## Residuals                703.35 805                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lm2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-7.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-8.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-9.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-10.png)<!-- -->

``` r
#Compare models
anova(lmer1, lm2) #Random effect improves model fit to data
```

    ## Data: legume_range_df_introducedY
    ## Models:
    ## lm2: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + scale(total_area_native) + annual + uses_num_uses
    ## lmer1: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + scale(total_area_native) + annual + uses_num_uses + (1 | genus)
    ##       npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)    
    ## lm2      9 2207.4 2249.7 -1094.7   2189.4                        
    ## lmer1   10 2192.8 2239.8 -1086.4   2172.8 16.66  1  4.472e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Introduced area (if introduced at all)
lmer2 <- lmer(log(total_area_introduced/1e+6)~EFN*Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+uses_num_uses+(1|genus), data=legume_range_df_introducedY) 
summary(lmer2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_introduced/1e+06) ~ EFN * Domatia + fixer + scale(abs_lat_native) +  
    ##     scale(total_area_native) + annual + uses_num_uses + (1 |      genus)
    ##    Data: legume_range_df_introducedY
    ## 
    ## REML criterion at convergence: 3701.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.0452 -0.4814  0.2038  0.6541  1.9387 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.8771   0.9366  
    ##  Residual             5.0312   2.2430  
    ## Number of obs: 813, groups:  genus, 195
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)              12.92192    0.34491  37.464
    ## EFN1                      0.25252    0.26108   0.967
    ## Domatia1                  1.46892    1.61220   0.911
    ## fixer1                   -0.37644    0.33906  -1.110
    ## scale(abs_lat_native)     0.36678    0.10371   3.537
    ## scale(total_area_native) -0.30987    0.09374  -3.306
    ## annual                    0.72724    0.27767   2.619
    ## uses_num_uses             0.52997    0.04187  12.657
    ## EFN1:Domatia1            -4.03588    2.09504  -1.926
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 scl(b__) scl(t__) annual uss_n_
    ## EFN1        -0.088                                                     
    ## Domatia1    -0.007  0.091                                              
    ## fixer1      -0.875  0.027  0.000                                       
    ## scl(bs_lt_)  0.104  0.106  0.034 -0.114                                
    ## scl(ttl_r_)  0.219 -0.031 -0.046 -0.129  0.033                         
    ## annual      -0.036  0.041 -0.005 -0.078 -0.078    0.132                
    ## uses_num_ss -0.384 -0.060 -0.004  0.098  0.044   -0.381    0.048       
    ## EFN1:Domat1  0.032 -0.138 -0.764 -0.021  0.009    0.064    0.002 -0.030

``` r
Anova(lmer2, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total_area_introduced/1e+06)
    ##                              Chisq Df Pr(>Chisq)    
    ## (Intercept)              1403.5566  1  < 2.2e-16 ***
    ## EFN                         0.9355  1  0.3334439    
    ## Domatia                     0.8302  1  0.3622270    
    ## fixer                       1.2327  1  0.2668898    
    ## scale(abs_lat_native)      12.5084  1  0.0004051 ***
    ## scale(total_area_native)   10.9277  1  0.0009474 ***
    ## annual                      6.8598  1  0.0088157 ** 
    ## uses_num_uses             160.2094  1  < 2.2e-16 ***
    ## EFN:Domatia                 3.7110  1  0.0540545 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-11.png)<!-- -->

``` r
#Native range size
lmer3 <- lmer(log(total_area_native/1e+6)~EFN*Domatia+fixer+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=df)
summary(lmer3)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(total_area_native/1e+06) ~ EFN * Domatia + fixer + annual +  
    ##     scale(abs_lat_native) + uses_num_uses + (1 | genus)
    ##    Data: df
    ## 
    ## REML criterion at convergence: 11273
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -9.6239 -0.3519  0.0649  0.5794  2.8799 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.3637   0.6031  
    ##  Residual             1.2914   1.1364  
    ## Number of obs: 3532, groups:  genus, 440
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           14.53436    0.09353 155.403
    ## EFN1                  -0.04090    0.08890  -0.460
    ## Domatia1               0.64120    0.32396   1.979
    ## fixer1                 0.09687    0.09716   0.997
    ## annual                 0.07016    0.06804   1.031
    ## scale(abs_lat_native) -0.12519    0.02789  -4.489
    ## uses_num_uses          0.28072    0.01408  19.943
    ## EFN1:Domatia1         -2.21933    0.67267  -3.299
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 annual sc(__) uss_n_
    ## EFN1        -0.018                                          
    ## Domatia1    -0.035  0.038                                   
    ## fixer1      -0.870 -0.007  0.019                            
    ## annual      -0.021 -0.007  0.008 -0.056                     
    ## scl(bs_lt_)  0.185  0.048  0.067 -0.167  0.062              
    ## uses_num_ss -0.166 -0.139  0.011  0.035  0.083  0.045       
    ## EFN1:Domat1  0.034 -0.111 -0.477 -0.026 -0.002  0.000 -0.031

``` r
Anova(lmer3, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total_area_native/1e+06)
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)           24150.0573  1  < 2.2e-16 ***
    ## EFN                       0.2116  1  0.6454956    
    ## Domatia                   3.9175  1  0.0477859 *  
    ## fixer                     0.9940  1  0.3187678    
    ## annual                    1.0632  1  0.3024768    
    ## scale(abs_lat_native)    20.1522  1  7.152e-06 ***
    ## uses_num_uses           397.7344  1  < 2.2e-16 ***
    ## EFN:Domatia              10.8853  1  0.0009693 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer3)
```

![](README_files/figure-gfm/Legume%20mixed%20models-12.png)<!-- -->

``` r
#Mycorrhizae
#Successful introduction?
binomial3 <- glmer(introducedY~AM*EM+scale(total_area_native)+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=subset(df, !is.na(myco)), family="binomial", control=glmerControl(optimizer="nloptwrap", optCtrl=list(maxfun=2e6)))
summary(binomial3) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ AM * EM + scale(total_area_native) + annual + scale(abs_lat_native) +  
    ##     uses_num_uses + (1 | genus)
    ##    Data: subset(df, !is.na(myco))
    ## Control: glmerControl(optimizer = "nloptwrap", optCtrl = list(maxfun = 2e+06))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    623.9    664.6   -303.0    605.9      672 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2264 -0.4656 -0.2142  0.4377  3.6810 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 0.8571   0.9258  
    ## Number of obs: 681, groups:  genus, 235
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -1.28855    0.51092  -2.522   0.0117 *  
    ## AMY                      -0.61326    0.50805  -1.207   0.2274    
    ## EMY                      -1.57862    0.75501  -2.091   0.0365 *  
    ## scale(total_area_native)  0.21653    0.14240   1.521   0.1284    
    ## annual                    0.46713    0.46861   0.997   0.3188    
    ## scale(abs_lat_native)     0.07854    0.13036   0.602   0.5469    
    ## uses_num_uses             0.98250    0.09709  10.119   <2e-16 ***
    ## AMY:EMY                   1.38056    0.82112   1.681   0.0927 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    scl(t__) annual scl(b__) uss_n_
    ## AMY         -0.898                                              
    ## EMY         -0.606  0.631                                       
    ## scl(ttl_r_)  0.076 -0.009  0.066                                
    ## annual      -0.146  0.074  0.087  0.085                         
    ## scl(bs_lt_) -0.040 -0.011 -0.020 -0.068   -0.061                
    ## uses_num_ss -0.224 -0.105 -0.088 -0.228    0.031  0.186         
    ## AMY:EMY      0.559 -0.620 -0.883 -0.031   -0.058  0.033    0.097

``` r
Anova(binomial3, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)                6.3606  1    0.01167 *  
    ## AM                         1.4570  1    0.22740    
    ## EM                         4.3717  1    0.03654 *  
    ## scale(total_area_native)   2.3121  1    0.12837    
    ## annual                     0.9937  1    0.31884    
    ## scale(abs_lat_native)      0.3630  1    0.54685    
    ## uses_num_uses            102.4036  1    < 2e-16 ***
    ## AM:EM                      2.8269  1    0.09270 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial3)
```

![](README_files/figure-gfm/Legume%20mixed%20models-13.png)<!-- -->

``` r
#Number of introduced ranges
lmer4 <- lmer(log(num_introduced)~AM+EM+scale(total_area_native)+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=subset(legume_range_df_introducedY, !is.na(myco)))
summary(lmer4)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(num_introduced) ~ AM + EM + scale(total_area_native) + annual +  
    ##     scale(abs_lat_native) + uses_num_uses + (1 | genus)
    ##    Data: subset(legume_range_df_introducedY, !is.na(myco))
    ## 
    ## REML criterion at convergence: 932.5
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.89248 -0.76275  0.08872  0.70696  2.06758 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.1001   0.3164  
    ##  Residual             0.9515   0.9754  
    ## Number of obs: 320, groups:  genus, 135
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.52812    0.23347   2.262
    ## AMY                       0.25360    0.22606   1.122
    ## EMY                      -0.02419    0.21793  -0.111
    ## scale(total_area_native) -0.13332    0.06181  -2.157
    ## annual                    0.42926    0.21861   1.964
    ## scale(abs_lat_native)    -0.04593    0.06272  -0.732
    ## uses_num_uses             0.31527    0.02571  12.260
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    scl(t__) annual scl(b__)
    ## AMY         -0.873                                       
    ## EMY         -0.181  0.117                                
    ## scl(ttl_r_)  0.058  0.017  0.095                         
    ## annual      -0.130  0.073  0.082  0.148                  
    ## scl(bs_lt_) -0.006 -0.012 -0.116 -0.062   -0.065         
    ## uses_num_ss -0.277 -0.119  0.006 -0.234   -0.026  0.105

``` r
Anova(lmer4, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)                5.1170  1    0.02369 *  
    ## AM                         1.2585  1    0.26194    
    ## EM                         0.0123  1    0.91163    
    ## scale(total_area_native)   4.6517  1    0.03102 *  
    ## annual                     3.8558  1    0.04957 *  
    ## scale(abs_lat_native)      0.5364  1    0.46393    
    ## uses_num_uses            150.3136  1    < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer4)
```

![](README_files/figure-gfm/Legume%20mixed%20models-14.png)<!-- -->

``` r
#Introduced area (if introduced at all)
lmer5 <- lmer(log(total_area_introduced/1e+6)~AM+EM+scale(total_area_native)+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=subset(legume_range_df_introducedY, !is.na(myco)))
summary(lmer5)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_introduced/1e+06) ~ AM + EM + scale(total_area_native) +  
    ##     annual + scale(abs_lat_native) + uses_num_uses + (1 | genus)
    ##    Data: subset(legume_range_df_introducedY, !is.na(myco))
    ## 
    ## REML criterion at convergence: 1467.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6881 -0.4034  0.2623  0.5761  1.8311 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 1.022    1.011   
    ##  Residual             4.952    2.225   
    ## Number of obs: 320, groups:  genus, 135
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)              11.95038    0.54748  21.828
    ## AMY                       0.97834    0.52523   1.863
    ## EMY                       0.09424    0.51997   0.181
    ## scale(total_area_native) -0.40625    0.14670  -2.769
    ## annual                    1.37971    0.51526   2.678
    ## scale(abs_lat_native)     0.54271    0.15244   3.560
    ## uses_num_uses             0.50075    0.06026   8.310
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    scl(t__) annual scl(b__)
    ## AMY         -0.865                                       
    ## EMY         -0.163  0.108                                
    ## scl(ttl_r_)  0.064  0.018  0.073                         
    ## annual      -0.128  0.077  0.064  0.170                  
    ## scl(bs_lt_) -0.008 -0.011 -0.115 -0.039   -0.027         
    ## uses_num_ss -0.284 -0.115 -0.001 -0.236   -0.028  0.103

``` r
Anova(lmer5, type=3, test="F")
```

    ## Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
    ## 
    ## Response: log(total_area_introduced/1e+06)
    ##                                 F Df Df.res    Pr(>F)    
    ## (Intercept)              472.8158  1 310.94 < 2.2e-16 ***
    ## AM                         3.4481  1 288.48 0.0643448 .  
    ## EM                         0.0319  1 309.32 0.8583800    
    ## scale(total_area_native)   7.5488  1 307.26 0.0063596 ** 
    ## annual                     7.0661  1 312.81 0.0082596 ** 
    ## scale(abs_lat_native)     12.4371  1 219.99 0.0005121 ***
    ## uses_num_uses             68.4783  1 309.62 3.867e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer5)
```

![](README_files/figure-gfm/Legume%20mixed%20models-15.png)<!-- -->

``` r
#Native range size
lmer6 <- lmer(log(total_area_native/1e+6)~AM+EM+annual+scale(abs_lat_native)+uses_num_uses+(1|genus), data=subset(df, !is.na(myco)))
summary(lmer6)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_native/1e+06) ~ AM + EM + annual + scale(abs_lat_native) +  
    ##     uses_num_uses + (1 | genus)
    ##    Data: subset(df, !is.na(myco))
    ## 
    ## REML criterion at convergence: 2068.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.0133 -0.3949  0.1228  0.6139  2.0675 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.2992   0.547   
    ##  Residual             1.0081   1.004   
    ## Number of obs: 681, groups:  genus, 235
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           15.13147    0.15637  96.765
    ## AMY                    0.09680    0.15020   0.644
    ## EMY                    0.04267    0.13685   0.312
    ## annual                -0.40833    0.17030  -2.398
    ## scale(abs_lat_native)  0.07115    0.05120   1.390
    ## uses_num_uses          0.15595    0.01951   7.995
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    annual sc(__)
    ## AMY         -0.873                            
    ## EMY         -0.362  0.267                     
    ## annual      -0.111  0.038  0.063              
    ## scl(bs_lt_)  0.007  0.000  0.007 -0.014       
    ## uses_num_ss -0.204 -0.073  0.056  0.020  0.109

``` r
Anova(lmer6, type=3, test="F")
```

    ## Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
    ## 
    ## Response: log(total_area_native/1e+06)
    ##                               F Df Df.res    Pr(>F)    
    ## (Intercept)           9322.9263  1 666.77 < 2.2e-16 ***
    ## AM                       0.4134  1 666.27   0.52044    
    ## EM                       0.0966  1 671.23   0.75610    
    ## annual                   5.7252  1 651.08   0.01701 *  
    ## scale(abs_lat_native)    1.9128  1 424.70   0.16738    
    ## uses_num_uses           63.6173  1 674.93 6.479e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer6)
```

![](README_files/figure-gfm/Legume%20mixed%20models-16.png)<!-- -->

#### Multiple mutualisms

We can also look at the effects of multiple mutualisms.

``` r
df$AM.num <- ifelse(df$AM == "Y", 1, ifelse(df$AM == "N", 0, NA))
df$EM.num <- ifelse(df$EM == "Y", 1, ifelse(df$EM == "N", 0, NA))

df$num_mutualisms <- as.numeric(as.character(df$fixer))+as.numeric(as.character(df$EFN))+as.numeric(as.character(df$AM.num))+as.numeric(as.character(df$Domatia))+as.numeric(as.character(df$EM.num))

summary.mnum2 <- ungroup(subset(df, !is.na(num_mutualisms)) %>% group_by(num_mutualisms, introducedY) %>% summarize(n=n()))
summary.mnum2.wide <- spread(summary.mnum2, key = introducedY, value=n)
colnames(summary.mnum2.wide) <- c("num_mutualisms","Not_introduced",  "Introduced")
summary.mnum2.wide$total <- summary.mnum2.wide$Not_introduced+summary.mnum2.wide$Introduced
summary.mnum2.wide$prop.introduced <- summary.mnum2.wide$Introduced/(summary.mnum2.wide$total)
prop.mnum <- paste0(summary.mnum2.wide$Introduced, "/", summary.mnum2.wide$total)

inset_p_mnum <- ggplot(data=subset(summary.mnum2.wide, !is.na(num_mutualisms)), aes(x=num_mutualisms, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+xlab("Mutualisms (no.)")+ylab("Introduced (prop.)")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 2, y = 0.55, label = "**")+geom_text(aes(x=num_mutualisms, y=0.05, label=prop.mnum), color="white")

summary.mnum <- subset(df, !is.na(num_mutualisms) & introducedY == 1) %>% group_by(num_mutualisms) %>% summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n))
summary.mnum$num_mutualisms <- as.factor(summary.mnum$num_mutualisms)

p_num_mutualisms <- ggplot(data=summary.mnum, aes(x=num_mutualisms, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num_mutualisms, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+theme_cowplot()+ylab("Introduced ranges (no.)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num_mutualisms, y= 3, label=n))+annotate("text", x=3, y=15.5, label="p = 0.06")
p_num_mutualisms
```

![](README_files/figure-gfm/Number%20of%20mutualisms-1.png)<!-- -->

``` r
lmer7 <-  lmer(log(num_introduced)~num_mutualisms+annual+scale(abs_lat_native)+uses_num_uses+scale(total_area_native)+(1|genus), data=subset(df, !is.na(num_mutualisms) & introducedY == 1))
summary(lmer7)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(num_introduced) ~ num_mutualisms + annual + scale(abs_lat_native) +  
    ##     uses_num_uses + scale(total_area_native) + (1 | genus)
    ##    Data: subset(df, !is.na(num_mutualisms) & introducedY == 1)
    ## 
    ## REML criterion at convergence: 931.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7349 -0.7400  0.1102  0.7065  2.0530 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.1419   0.3767  
    ##  Residual             0.9134   0.9557  
    ## Number of obs: 320, groups:  genus, 135
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.38490    0.22409   1.718
    ## num_mutualisms            0.18122    0.09737   1.861
    ## annual                    0.39816    0.21709   1.834
    ## scale(abs_lat_native)    -0.05742    0.06321  -0.908
    ## uses_num_uses             0.31697    0.02540  12.478
    ## scale(total_area_native) -0.13926    0.06170  -2.257
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_mtl annual scl(b__) uss_n_
    ## num_mutlsms -0.862                              
    ## annual      -0.087  0.030                       
    ## scl(bs_lt_)  0.008 -0.040 -0.037                
    ## uses_num_ss -0.358 -0.045 -0.021  0.107         
    ## scl(ttl_r_)  0.114 -0.033  0.154 -0.039   -0.235

``` r
Anova(lmer7, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)                2.9504  1    0.08586 .  
    ## num_mutualisms             3.4640  1    0.06272 .  
    ## annual                     3.3636  1    0.06665 .  
    ## scale(abs_lat_native)      0.8253  1    0.36364    
    ## uses_num_uses            155.7070  1    < 2e-16 ***
    ## scale(total_area_native)   5.0947  1    0.02400 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer7)
```

![](README_files/figure-gfm/Number%20of%20mutualisms-2.png)<!-- -->

``` r
binomial4 <-  glmer(introducedY~num_mutualisms+annual+scale(abs_lat_native)+uses_num_uses+scale(total_area_native)+(1|genus), data=subset(df, !is.na(num_mutualisms)), family="binomial")
summary(binomial4)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: introducedY ~ num_mutualisms + annual + scale(abs_lat_native) +  
    ##     uses_num_uses + scale(total_area_native) + (1 | genus)
    ##    Data: subset(df, !is.na(num_mutualisms))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    617.9    649.6   -302.0    603.9      674 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5352 -0.4407 -0.2400  0.4115  3.4945 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 0.8366   0.9147  
    ## Number of obs: 681, groups:  genus, 235
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.93749    0.46314  -6.343 2.26e-10 ***
    ## num_mutualisms            0.51364    0.19944   2.575    0.010 *  
    ## annual                    0.59541    0.46719   1.274    0.203    
    ## scale(abs_lat_native)     0.04971    0.13036   0.381    0.703    
    ## uses_num_uses             0.98389    0.09650  10.196  < 2e-16 ***
    ## scale(total_area_native)  0.22095    0.14284   1.547    0.122    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_mtl annual scl(b__) uss_n_
    ## num_mutlsms -0.884                              
    ## annual      -0.105  0.029                       
    ## scl(bs_lt_)  0.005 -0.058 -0.073                
    ## uses_num_ss -0.439  0.107  0.042  0.165         
    ## scl(ttl_r_)  0.144 -0.057  0.063 -0.065   -0.252

``` r
Anova(binomial4, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               40.2283  1  2.259e-10 ***
    ## num_mutualisms             6.6329  1    0.01001 *  
    ## annual                     1.6242  1    0.20251    
    ## scale(abs_lat_native)      0.1454  1    0.70294    
    ## uses_num_uses            103.9575  1  < 2.2e-16 ***
    ## scale(total_area_native)   2.3925  1    0.12191    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial4)
```

![](README_files/figure-gfm/Number%20of%20mutualisms-3.png)<!-- -->

``` r
fig3 <- plot_grid(inset_p_mnum, p_num_mutualisms, ncol=1, labels="AUTO")
fig3
```

![](README_files/figure-gfm/Number%20of%20mutualisms-4.png)<!-- -->

``` r
save_plot("Figure3.pdf", fig3, base_height = 8, base_width=5)

sum.nmum <- subset(df, !is.na(num_mutualisms)) %>% group_by(AM, EM, Domatia, EFN, fixer, num_mutualisms) %>% summarize(n=n())
kable(sum.nmum)
```

| AM  | EM  | Domatia | EFN | fixer | num_mutualisms |   n |
|:----|:----|:--------|:----|:------|---------------:|----:|
| N   | N   | 0       | 0   | 0     |              0 |   4 |
| N   | N   | 0       | 0   | 1     |              1 |  28 |
| N   | N   | 0       | 1   | 0     |              1 |   1 |
| N   | Y   | 0       | 0   | 0     |              1 |  19 |
| N   | Y   | 0       | 0   | 1     |              2 |  16 |
| N   | Y   | 0       | 1   | 1     |              3 |   3 |
| Y   | N   | 0       | 0   | 0     |              1 |  63 |
| Y   | N   | 0       | 0   | 1     |              2 | 430 |
| Y   | N   | 0       | 1   | 0     |              2 |  13 |
| Y   | N   | 0       | 1   | 1     |              3 |  63 |
| Y   | N   | 1       | 0   | 0     |              2 |   1 |
| Y   | N   | 1       | 0   | 1     |              3 |   2 |
| Y   | Y   | 0       | 0   | 0     |              2 |  13 |
| Y   | Y   | 0       | 0   | 1     |              3 |  52 |
| Y   | Y   | 0       | 1   | 1     |              4 |  14 |
| Y   | Y   | 1       | 0   | 0     |              3 |   1 |

## Ant dataset

We obtained introduced and native range size data for ants from Benoit
Guenard and Evan Economo, and trait data on which ant species visit
EFNs, disperse seeds, and nest in domatia from Kaur et al. (2019).

``` r
###Native range size
antarea <- read.csv("Ant_species_native range.csv") #native area for ants
invarea <- read.csv("invaded_area_ant_species.csv") #invaded area for ants - all types of introduction and all available data
abslat <- read.csv("absolute_native_lat_ants7Feb.csv") #Absolute midpoint latitude of native range

#merging invaded and native area datasets
nat_inv_area <- merge(invarea, antarea, by='Phy', all.y=TRUE) 
nat_inv_area$total.area.introduced <- ifelse(is.na(nat_inv_area$total.area.introduced), 0, nat_inv_area$total.area.introduced) #Make NAs zeros 
nat_inv_area_lat <- merge(abslat[, c("Phy", "abs_lat_native")], nat_inv_area, by = 'Phy') #Add midpoint native latitude

#Trait data from Kaur et al. (2019)
antefn <- read.csv("Species_EFN_Data.csv")
antdom <- read.csv("Species_Domatia_Data.csv")
antseed <- read.csv("Species_Seed_Dispersal_Data.csv")

#Creating merged datasets for area
efn_area <- merge(nat_inv_area_lat, antefn, by.y = "Phy", all = TRUE)
efn_dom_area <- merge(efn_area, antdom, by.y = "Phy", all = TRUE)
area <- merge(efn_dom_area, antseed, by.y = "Phy", all = TRUE)

#Make sure factors are factors
area$EFN <- as.factor(area$EFN)
area$Seed_Dispersal <- as.factor(area$Seed_Dispersal)
area$Domatia <- as.factor(area$Domatia)

#Exotic, intercepted, or indoor introduced
exotic <- read.csv("alien_exotic.csv") 
indoor <- read.csv("alien_indoor.csv")
intercepted <- read.csv("alien_intercepted.csv")

#Add exotic, intercepted, or indoor introduced
area <- merge(area, exotic[, c("Phy", "Exotic_status")], by="Phy", all=TRUE)
colnames(area)[[14]] <- "ExoticY"
area <- merge(area, indoor[, c("Phy", "Exotic_status")], by="Phy", all=TRUE)
colnames(area)[[15]] <- "IndoorY"
area <- merge(area, intercepted[, c("Phy", "Exotic_status")], by="Phy", all=TRUE)
colnames(area)[[16]] <- "InterceptedY"

#Introduced? 
area$introducedY <- ifelse(area$total.area.introduced > 0, 1, 0)
area$IndoorY <- ifelse(!is.na(area$IndoorY) & area$IndoorY == "Indoor Introduced" & area$introducedY == 1, 1, ifelse(is.na(area$IndoorY) & area$introducedY == 1, 0, NA))
area$ExoticY <- ifelse(!is.na(area$ExoticY) & area$ExoticY == "Exotic" & area$introducedY == 1, 1, ifelse(is.na(area$ExoticY) & area$introducedY == 1, 0, NA))
area$InterceptedY <- ifelse(!is.na(area$InterceptedY) & area$InterceptedY == "Intercepted" & area$introducedY == 1, 1, ifelse(is.na(area$InterceptedY) & area$introducedY == 1, 0, NA))
```

### Summarize ant dataset

``` r
sum(complete.cases(subset(area, !is.na(introducedY))[, c("EFN", "Domatia", "Seed_Dispersal")])) #Total taxa with at least some trait data 
```

    ## [1] 3023

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(EFN) %>% summarize(n=n()))
```

| EFN |     n |
|:----|------:|
| 0   |  2910 |
| 1   |   115 |
| NA  | 12271 |

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(Domatia) %>% summarize(n=n()))
```

| Domatia |     n |
|:--------|------:|
| 0       |  2967 |
| 1       |    58 |
| NA      | 12271 |

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(Seed_Dispersal) %>% summarize(n=n()))
```

| Seed_Dispersal |     n |
|:---------------|------:|
| 0              |  2758 |
| 1              |   297 |
| NA             | 12241 |

``` r
pt_size <- 3
y_limits <- c(-100000, 1.5e+7)
er_width <- 0.1
y_text <- -100000

summary.ant.efn <- subset(area, introducedY == 1 & !is.na(EFN)) %>% group_by(EFN) %>% summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antEFN <- ggplot(data=summary.ant.efn, aes(x=EFN, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range (sq. km)")+xlab("Visits EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "*")

summary.ant.efn2 <- ungroup(subset(area, !is.na(EFN)) %>% group_by(EFN, introducedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'EFN'. You can override using the `.groups`
    ## argument.

``` r
summary.ant.efn2.wide <- spread(summary.ant.efn2, key = introducedY, value=n)
colnames(summary.ant.efn2.wide) <- c("EFN","Not_introduced",  "Introduced")
summary.ant.efn2.wide$total <- summary.ant.efn2.wide$Not_introduced+summary.ant.efn2.wide$Introduced
summary.ant.efn2.wide$prop.introduced <- summary.ant.efn2.wide$Introduced/(summary.ant.efn2.wide$total)
prop.ant.efn <- paste0(summary.ant.efn2.wide$Introduced, "/", summary.ant.efn2.wide$total)

inset_p_antEFN <- ggplot(data=summary.ant.efn2.wide, aes(x=EFN, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=EFN, y=0.05, label=prop.ant.efn), color="white")

summary.ant.dom <- subset(area, introducedY == 1 & !is.na(Domatia)) %>% group_by(Domatia) %>% summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antdom <- ggplot(data=summary.ant.dom, aes(x=Domatia, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range (sq. km)")+xlab("Nests in domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "ns")

summary.ant.dom2 <- ungroup(subset(area, !is.na(Domatia)) %>% group_by(Domatia, introducedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Domatia'. You can override using the
    ## `.groups` argument.

``` r
summary.ant.dom2.wide <- spread(summary.ant.dom2, key = introducedY, value=n)
colnames(summary.ant.dom2.wide) <- c("Domatia","Not_introduced",  "Introduced")
summary.ant.dom2.wide$total <- summary.ant.dom2.wide$Not_introduced+summary.ant.dom2.wide$Introduced
summary.ant.dom2.wide$prop.introduced <- summary.ant.dom2.wide$Introduced/(summary.ant.dom2.wide$total)
prop.ant.dom <- paste0(summary.ant.dom2.wide$Introduced, "/", summary.ant.dom2.wide$total)

inset_p_antdom <- ggplot(data=summary.ant.dom2.wide, aes(x=Domatia, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.ant.dom), color="white")

summary.ant.elaiosome <- subset(area, introducedY == 1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal) %>% summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_elaiosome <- ggplot(data=summary.ant.elaiosome, aes(x=Seed_Dispersal, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Seed_Dispersal, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range (sq. km)")+xlab("Disperses seeds")+geom_text(aes(x=Seed_Dispersal, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "ns")

summary.ant.seed2 <- ungroup(subset(area, !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, introducedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Seed_Dispersal'. You can override using
    ## the `.groups` argument.

``` r
summary.ant.seed2.wide <- spread(summary.ant.seed2, key = introducedY, value=n)
colnames(summary.ant.seed2.wide) <- c("Seed_Dispersal","Not_introduced",  "Introduced")
summary.ant.seed2.wide$total <- summary.ant.seed2.wide$Not_introduced+summary.ant.seed2.wide$Introduced
summary.ant.seed2.wide$prop.introduced <- summary.ant.seed2.wide$Introduced/(summary.ant.seed2.wide$total)
prop.ant.seed <- paste0(summary.ant.seed2.wide$Introduced, "/", summary.ant.seed2.wide$total)

inset_p_antseed <- ggplot(data=summary.ant.seed2.wide, aes(x=Seed_Dispersal, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.ant.seed), color="white")

fig4top <- plot_grid(inset_p_antEFN, inset_p_antseed, inset_p_antdom,  nrow=1, labels="AUTO")
fig4bottom <- plot_grid(p_antEFN, p_elaiosome, p_antdom, nrow=1, labels=c("D", "E", "F"))
fig4 <- plot_grid(fig4top, fig4bottom, nrow=2)
fig4
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
save_plot("Figure4.pdf", fig4, base_width=8, base_height =8)
```

#### Native range size

``` r
y_limits = c(3200000, 1.25e+07)
y_text = 3300000
summary.ant.efn.native <- subset(area, !is.na(EFN) & !is.na(total.area.native)) %>% group_by(EFN) %>% summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antEFN_native <- ggplot(data=summary.ant.efn.native, aes(x=EFN, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Visits EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 8100000, label = "***")

summary.ant.dom.native <- subset(area, !is.na(Domatia)& !is.na(total.area.native)) %>% group_by(Domatia) %>% summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antdom_native <- ggplot(data=summary.ant.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nests in domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 8100000, label = "*")

summary.ant.seed.native <- subset(area, !is.na(Seed_Dispersal)& !is.na(total.area.native)) %>% group_by(Seed_Dispersal) %>% summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antseed_native <- ggplot(data=summary.ant.seed.native, aes(x=Seed_Dispersal, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Seed_Dispersal, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Disperses seeds")+geom_text(aes(x=Seed_Dispersal, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 8100000, label = "***")

fig5 <- plot_grid(p_antEFN_native, p_antseed_native, p_antdom_native, nrow=1, labels="AUTO")
fig5
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
save_plot("Figure5.pdf", fig5, base_height = 4, base_width = 12)
```

### Multiple mutualisms

``` r
area$genus <- as.factor(word(area$Phy, 1, 1, sep="_")) #Extract ant genus
area$num.mm <- as.numeric(as.character(area$Seed_Dispersal))+as.numeric(as.character(area$Domatia))+as.numeric(as.character(area$EFN))

binomial10 <- glmer(introducedY ~ num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=area, family ="binomial")
summary(binomial10)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ num.mm + scale(abs_lat_native) + scale(total.area.native) +  
    ##     (1 | genus)
    ##    Data: area
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1742.6   1772.6   -866.3   1732.6     2991 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.4583 -0.3226 -0.2144 -0.1310 15.9280 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 1.604    1.267   
    ## Number of obs: 2996, groups:  genus, 249
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -3.43863    0.20858 -16.486  < 2e-16 ***
    ## num.mm                    0.95019    0.14016   6.780 1.21e-11 ***
    ## scale(abs_lat_native)    -0.17249    0.07052  -2.446   0.0144 *  
    ## scale(total.area.native)  0.45390    0.04010  11.319  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) num.mm sc(__)
    ## num.mm      -0.155              
    ## scl(bs_lt_) -0.021  0.004       
    ## scl(ttl.r.) -0.322 -0.203 -0.079

``` r
Anova(binomial10, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)              271.7739  1  < 2.2e-16 ***
    ## num.mm                    45.9628  1  1.205e-11 ***
    ## scale(abs_lat_native)      5.9824  1    0.01445 *  
    ## scale(total.area.native) 128.1307  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial10)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
lmer11 <- lmer(log(total.area.introduced)~num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=subset(area, introducedY == 1))
summary(lmer11)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(total.area.introduced) ~ num.mm + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | genus)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ## REML criterion at convergence: 1528.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.9381 -0.6164  0.0664  0.6422  2.0419 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.8491   0.9214  
    ##  Residual             4.6031   2.1455  
    ## Number of obs: 341, groups:  genus, 69
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               12.4112     0.2018  61.514
    ## num.mm                     0.3827     0.2034   1.882
    ## scale(abs_lat_native)     -0.3843     0.1276  -3.012
    ## scale(total.area.native)   0.2076     0.1218   1.705
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) num.mm sc(__)
    ## num.mm      -0.404              
    ## scl(bs_lt_)  0.013 -0.051       
    ## scl(ttl.r.)  0.055 -0.359 -0.054

``` r
Anova(lmer11, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total.area.introduced)
    ##                              Chisq Df Pr(>Chisq)    
    ## (Intercept)              3783.9495  1  < 2.2e-16 ***
    ## num.mm                      3.5413  1   0.059859 .  
    ## scale(abs_lat_native)       9.0728  1   0.002594 ** 
    ## scale(total.area.native)    2.9064  1   0.088229 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary.ant.mm <- ungroup(subset(area, !is.na(introducedY)) %>% group_by(num.mm, introducedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'num.mm'. You can override using the
    ## `.groups` argument.

``` r
summary.ant.mm.wide <- spread(summary.ant.mm, key=introducedY, value=n)
colnames(summary.ant.mm.wide) <- c("Num_mutualisms","Not_introduced",  "Introduced")
summary.ant.mm.wide$total <- summary.ant.mm.wide$Not_introduced+summary.ant.mm.wide$Introduced
summary.ant.mm.wide$prop.introduced <- summary.ant.mm.wide$Introduced/(summary.ant.mm.wide$total)
prop.ant.mm <- paste0(summary.ant.mm.wide$Introduced, "/", summary.ant.mm.wide$total)
prop.ant.mm <- prop.ant.mm[1:4]

inset_p_antmm <- ggplot(data=subset(summary.ant.mm.wide, !is.na(Num_mutualisms)),aes(x=Num_mutualisms, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+ylab("Introduced (prop.)")+xlab("Mutualisms (no.)")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=Num_mutualisms, y=0.05, label=prop.ant.mm), color="white")

summary.ant.mnum <- subset(area, !is.na(num.mm) & introducedY == 1) %>% group_by(num.mm) %>% summarize(n=n(), mean_area_introduced = mean(total.area.introduced, na.rm=TRUE), sd_area_introduced = sd(total.area.introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n))
summary.ant.mnum$num.mm <- as.factor(summary.ant.mnum$num.mm)

p_num_mm <- ggplot(data=summary.ant.mnum, aes(x=num.mm, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num.mm, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+theme_cowplot()+ylab("Introduced range (sq. km)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num.mm, y= -1000000, label=n))+scale_y_continuous(limits=c(-3000000, 3.5e+7))+annotate("text", x=2.5, y=2e+7, label="p = 0.06")

fig6 <- plot_grid(inset_p_antmm, p_num_mm, ncol=1, labels="AUTO")
fig6
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
save_plot("Figure6.pdf", fig6, base_height = 8, base_width=5)
```

### Mixed models

``` r
#Introduction success
binomial6 <- glmer(introducedY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=area, family="binomial")
summary(binomial6)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | genus)
    ##    Data: area
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1733.8   1775.8   -859.9   1719.8     2989 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.2146 -0.3179 -0.2098 -0.1300 16.6967 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 1.552    1.246   
    ## Number of obs: 2996, groups:  genus, 249
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               -3.4136     0.2062 -16.553  < 2e-16 ***
    ## EFN1                       0.9321     0.2689   3.467 0.000527 ***
    ## Domatia1                  -0.5262     0.5147  -1.022 0.306683    
    ## Seed_Dispersal1            1.2882     0.2026   6.357 2.05e-10 ***
    ## scale(abs_lat_native)     -0.2142     0.0719  -2.978 0.002897 ** 
    ## scale(total.area.native)   0.4451     0.0404  11.018  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 Sd_Ds1 sc(__)
    ## EFN1        -0.052                            
    ## Domatia1    -0.037 -0.115                     
    ## Sed_Dsprsl1 -0.127 -0.114 -0.004              
    ## scl(bs_lt_) -0.009  0.064  0.095 -0.135       
    ## scl(ttl.r.) -0.326 -0.076 -0.051 -0.175 -0.068

``` r
Anova(binomial6, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)              273.9997  1  < 2.2e-16 ***
    ## EFN                       12.0170  1  0.0005272 ***
    ## Domatia                    1.0449  1  0.3066831    
    ## Seed_Dispersal            40.4152  1  2.053e-10 ***
    ## scale(abs_lat_native)      8.8712  1  0.0028971 ** 
    ## scale(total.area.native) 121.3918  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial6)
```

![](README_files/figure-gfm/ant%20glmms-1.png)<!-- -->

``` r
#Total introduced area
lmer9 <- lmer(log(total.area.introduced)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=subset(area, introducedY == 1))
summary(lmer9)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(total.area.introduced) ~ EFN + Domatia + Seed_Dispersal +  
    ##     scale(abs_lat_native) + scale(total.area.native) + (1 | genus)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ## REML criterion at convergence: 1524.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.9631 -0.6036  0.0588  0.5978  2.0171 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.819    0.905   
    ##  Residual             4.616    2.148   
    ## Number of obs: 341, groups:  genus, 69
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)              12.44926    0.20237  61.516
    ## EFN1                      0.83092    0.39157   2.122
    ## Domatia1                  0.25332    0.88489   0.286
    ## Seed_Dispersal1           0.05104    0.31413   0.162
    ## scale(abs_lat_native)    -0.34694    0.13149  -2.639
    ## scale(total.area.native)  0.22194    0.12247   1.812
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 Sd_Ds1 sc(__)
    ## EFN1        -0.134                            
    ## Domatia1    -0.035 -0.238                     
    ## Sed_Dsprsl1 -0.362 -0.200  0.032              
    ## scl(bs_lt_)  0.046  0.096  0.106 -0.208       
    ## scl(ttl.r.)  0.069 -0.135 -0.034 -0.305 -0.028

``` r
Anova(lmer9, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total.area.introduced)
    ##                              Chisq Df Pr(>Chisq)    
    ## (Intercept)              3784.2427  1  < 2.2e-16 ***
    ## EFN                         4.5029  1   0.033837 *  
    ## Domatia                     0.0820  1   0.774665    
    ## Seed_Dispersal              0.0264  1   0.870929    
    ## scale(abs_lat_native)       6.9621  1   0.008325 ** 
    ## scale(total.area.native)    3.2840  1   0.069959 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer9)
```

![](README_files/figure-gfm/ant%20glmms-2.png)<!-- -->

``` r
#Native range area
lmer10 <- lmer(log(total.area.native)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+(1|genus), data=area)
summary(lmer10)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     (1 | genus)
    ##    Data: area
    ## 
    ## REML criterion at convergence: 11466.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -6.6863 -0.4933  0.1685  0.6563  2.7779 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.522    0.7225  
    ##  Residual             2.506    1.5830  
    ## Number of obs: 2996, groups:  genus, 249
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           14.18328    0.06826 207.772
    ## EFN1                   0.56619    0.16131   3.510
    ## Domatia1               0.49247    0.22954   2.145
    ## Seed_Dispersal1        1.25268    0.11176  11.209
    ## scale(abs_lat_native)  0.06977    0.03362   2.075
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 Sd_Ds1
    ## EFN1        -0.052                     
    ## Domatia1    -0.073 -0.109              
    ## Sed_Dsprsl1 -0.099 -0.186 -0.003       
    ## scl(bs_lt_) -0.062  0.051  0.080 -0.073

``` r
Anova(lmer10, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total.area.native)
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)           43169.1453  1  < 2.2e-16 ***
    ## EFN                      12.3193  1  0.0004483 ***
    ## Domatia                   4.6030  1  0.0319168 *  
    ## Seed_Dispersal          125.6414  1  < 2.2e-16 ***
    ## scale(abs_lat_native)     4.3058  1  0.0379831 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer10)
```

![](README_files/figure-gfm/ant%20glmms-3.png)<!-- -->

### By introduction mode

``` r
# Exotic
area$ExoticY <- as.factor(area$ExoticY)
exotic.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, ExoticY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'EFN'. You can override using the `.groups`
    ## argument.

``` r
exotic.efn.wide <- spread(exotic.efn, key = ExoticY, value=n)
colnames(exotic.efn.wide) <- c("EFN","Not_exotic",  "Exotic")
exotic.efn.wide$EFN <- c("No", "Yes")
exotic.efn.wide$total <- exotic.efn.wide$Not_exotic + exotic.efn.wide$Exotic
exotic.efn.wide$prop.exotic <- exotic.efn.wide$Exotic/exotic.efn.wide$total
prop.exotic.efn <- paste0(exotic.efn.wide$Exotic, "/", exotic.efn.wide$total)

p_exoticefn <- ggplot(data=exotic.efn.wide, aes(x=EFN, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "p = 0.08")+geom_text(aes(x=EFN, y=0.05, label=prop.exotic.efn), color="white")

exotic.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, ExoticY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Domatia'. You can override using the
    ## `.groups` argument.

``` r
exotic.dom.wide <- spread(exotic.dom, key = ExoticY, value=n)
colnames(exotic.dom.wide) <- c("Domatia","Not_exotic",  "Exotic")
exotic.dom.wide$Domatia <- c("No", "Yes")
exotic.dom.wide$total <- exotic.dom.wide$Not_exotic + exotic.dom.wide$Exotic
exotic.dom.wide$prop.exotic <- exotic.dom.wide$Exotic/exotic.dom.wide$total
prop.exotic.dom <- paste0(exotic.dom.wide$Exotic, "/", exotic.dom.wide$total)

p_exoticdom <- ggplot(data=exotic.dom.wide, aes(x=Domatia, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.exotic.dom), color="white")

exotic.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, ExoticY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Seed_Dispersal'. You can override using
    ## the `.groups` argument.

``` r
exotic.seed.wide <- spread(exotic.seed, key = ExoticY, value=n)
colnames(exotic.seed.wide) <- c("Seed_Dispersal","Not_exotic",  "Exotic")
exotic.seed.wide$Seed_Dispersal <- c("No", "Yes")
exotic.seed.wide$total <- exotic.seed.wide$Not_exotic + exotic.seed.wide$Exotic
exotic.seed.wide$prop.exotic <- exotic.seed.wide$Exotic/exotic.seed.wide$total
prop.exotic.seed <- paste0(exotic.seed.wide$Exotic, "/", exotic.seed.wide$total)

p_exoticseed <- ggplot(data=exotic.seed.wide, aes(x=Seed_Dispersal, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.exotic.seed), color="white")

p_exotic <- plot_grid(p_exoticefn, p_exoticdom, p_exoticseed, nrow=1, labels=c("G", "H", "I"))

binomial7 <- glmer(ExoticY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=subset(area, introducedY == 1), family="binomial")
summary(binomial7)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: ExoticY ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | genus)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    419.3    446.1   -202.6    405.3      334 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.4743 -0.7732  0.3636  0.5996  1.7443 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 1.794    1.339   
    ## Number of obs: 341, groups:  genus, 69
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)               0.72411    0.26528   2.730  0.00634 **
    ## EFN1                      0.78177    0.43911   1.780  0.07502 . 
    ## Domatia1                 -0.41382    0.94543  -0.438  0.66160   
    ## Seed_Dispersal1          -0.19170    0.35916  -0.534  0.59351   
    ## scale(abs_lat_native)    -0.44570    0.16192  -2.753  0.00591 **
    ## scale(total.area.native) -0.01343    0.13406  -0.100  0.92020   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 Sd_Ds1 sc(__)
    ## EFN1        -0.072                            
    ## Domatia1    -0.028 -0.210                     
    ## Sed_Dsprsl1 -0.331 -0.215  0.021              
    ## scl(bs_lt_) -0.028  0.028  0.115 -0.171       
    ## scl(ttl.r.)  0.078 -0.116 -0.037 -0.324 -0.038

``` r
Anova(binomial7, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: ExoticY
    ##                           Chisq Df Pr(>Chisq)   
    ## (Intercept)              7.4504  1   0.006342 **
    ## EFN                      3.1697  1   0.075018 . 
    ## Domatia                  0.1916  1   0.661602   
    ## Seed_Dispersal           0.2849  1   0.593514   
    ## scale(abs_lat_native)    7.5770  1   0.005912 **
    ## scale(total.area.native) 0.0100  1   0.920201   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial7)
```

![](README_files/figure-gfm/Introduction%20mode-1.png)<!-- -->

``` r
#Indoors?
area$IndoorY <- as.factor(area$IndoorY)

indoor.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, IndoorY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'EFN'. You can override using the `.groups`
    ## argument.

``` r
indoor.efn.wide <- spread(indoor.efn, key = IndoorY, value=n)
colnames(indoor.efn.wide) <- c("EFN","Not_indoor",  "indoor")
indoor.efn.wide$EFN <- c("No", "Yes")
indoor.efn.wide$total <- indoor.efn.wide$Not_indoor + indoor.efn.wide$indoor
indoor.efn.wide$prop.indoor <- indoor.efn.wide$indoor/indoor.efn.wide$total
prop.indoor.efn <- paste0(indoor.efn.wide$indoor, "/", indoor.efn.wide$total)

p_indoorefn <- ggplot(data=indoor.efn.wide, aes(x=EFN, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "*")+geom_text(aes(x=EFN, y=0.05, label=prop.indoor.efn), color="white")

indoor.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, IndoorY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Domatia'. You can override using the
    ## `.groups` argument.

``` r
indoor.dom.wide <- spread(indoor.dom, key = IndoorY, value=n)
colnames(indoor.dom.wide) <- c("Domatia","Not_indoor",  "indoor")
indoor.dom.wide$Domatia <- c("No", "Yes")
indoor.dom.wide$total <- indoor.dom.wide$Not_indoor + indoor.dom.wide$indoor
indoor.dom.wide$prop.indoor <- indoor.dom.wide$indoor/indoor.dom.wide$total
prop.indoor.dom <- paste0(indoor.dom.wide$indoor, "/", indoor.dom.wide$total)

p_indoordom <- ggplot(data=indoor.dom.wide, aes(x=Domatia, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.indoor.dom), color="white")

indoor.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, IndoorY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Seed_Dispersal'. You can override using
    ## the `.groups` argument.

``` r
indoor.seed.wide <- spread(indoor.seed, key = IndoorY, value=n)
colnames(indoor.seed.wide) <- c("Seed_Dispersal","Not_indoor",  "indoor")
indoor.seed.wide$Seed_Dispersal <- c("No", "Yes")
indoor.seed.wide$total <- indoor.seed.wide$Not_indoor + indoor.seed.wide$indoor
indoor.seed.wide$prop.indoor <- indoor.seed.wide$indoor/indoor.seed.wide$total
prop.indoor.seed <- paste0(indoor.seed.wide$indoor, "/", indoor.seed.wide$total)

p_indoorseed <- ggplot(data=indoor.seed.wide, aes(x=Seed_Dispersal, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "**")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.indoor.seed), color="white")

p_indoor <- plot_grid(p_indoorefn, p_indoordom, p_indoorseed, nrow=1, labels=c("D","E", "F"))

binomial8 <- glmer(IndoorY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=subset(area, introducedY == 1), family="binomial")
summary(binomial8)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: IndoorY ~ EFN + Seed_Dispersal + Domatia + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | genus)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    384.5    411.3   -185.2    370.5      334 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7605 -0.5556 -0.3933  0.7043  3.0577 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 1.077    1.038   
    ## Number of obs: 341, groups:  genus, 69
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               -1.1898     0.2635  -4.516  6.3e-06 ***
    ## EFN1                       1.0667     0.4293   2.485  0.01296 *  
    ## Seed_Dispersal1           -1.0376     0.4007  -2.589  0.00962 ** 
    ## Domatia1                   0.2079     0.9904   0.210  0.83373    
    ## scale(abs_lat_native)     -0.1590     0.1523  -1.044  0.29650    
    ## scale(total.area.native)   0.3360     0.1378   2.438  0.01477 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Sd_Ds1 Domat1 sc(__)
    ## EFN1        -0.124                            
    ## Sed_Dsprsl1 -0.238 -0.277                     
    ## Domatia1    -0.074 -0.263  0.022              
    ## scl(bs_lt_)  0.060  0.079 -0.171  0.099       
    ## scl(ttl.r.)  0.003 -0.085 -0.362 -0.013 -0.048

``` r
Anova(binomial8, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: IndoorY
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)              20.3954  1  6.298e-06 ***
    ## EFN                       6.1745  1    0.01296 *  
    ## Seed_Dispersal            6.7039  1    0.00962 ** 
    ## Domatia                   0.0441  1    0.83373    
    ## scale(abs_lat_native)     1.0899  1    0.29650    
    ## scale(total.area.native)  5.9432  1    0.01477 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial8)
```

![](README_files/figure-gfm/Introduction%20mode-2.png)<!-- -->

``` r
#Intercepted
area$InterceptedY <- as.factor(area$InterceptedY)

intercepted.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, InterceptedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'EFN'. You can override using the `.groups`
    ## argument.

``` r
intercepted.efn.wide <- spread(intercepted.efn, key = InterceptedY, value=n)
colnames(intercepted.efn.wide) <- c("EFN","Not_intercepted",  "intercepted")
intercepted.efn.wide$EFN <- c("No", "Yes")
intercepted.efn.wide$total <- intercepted.efn.wide$Not_intercepted + intercepted.efn.wide$intercepted
intercepted.efn.wide$prop.intercepted <- intercepted.efn.wide$intercepted/intercepted.efn.wide$total
prop.intercepted.efn <- paste0(intercepted.efn.wide$intercepted, "/", intercepted.efn.wide$total)

p_interceptedefn <- ggplot(data=intercepted.efn.wide, aes(x=EFN, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=EFN, y=0.05, label=prop.intercepted.efn), color="white")

intercepted.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, InterceptedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Domatia'. You can override using the
    ## `.groups` argument.

``` r
intercepted.dom.wide <- spread(intercepted.dom, key = InterceptedY, value=n)
colnames(intercepted.dom.wide) <- c("Domatia","Not_intercepted",  "intercepted")
intercepted.dom.wide$Not_intercepted <- ifelse(is.na(intercepted.dom.wide$Not_intercepted), 0, intercepted.dom.wide$Not_intercepted)
intercepted.dom.wide$Domatia <- c("No", "Yes")
intercepted.dom.wide$total <- intercepted.dom.wide$Not_intercepted + intercepted.dom.wide$intercepted
intercepted.dom.wide$prop.intercepted <- intercepted.dom.wide$intercepted/intercepted.dom.wide$total
prop.intercepted.dom <- paste0(intercepted.dom.wide$intercepted, "/", intercepted.dom.wide$total)

p_intercepteddom <- ggplot(data=intercepted.dom.wide, aes(x=Domatia, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.45, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.intercepted.dom), color="white")

intercepted.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, InterceptedY) %>% summarize(n=n()))
```

    ## `summarise()` has grouped output by 'Seed_Dispersal'. You can override using
    ## the `.groups` argument.

``` r
intercepted.seed.wide <- spread(intercepted.seed, key = InterceptedY, value=n)
colnames(intercepted.seed.wide) <- c("Seed_Dispersal","Not_intercepted",  "intercepted")
intercepted.seed.wide$Seed_Dispersal <- c("No", "Yes")
intercepted.seed.wide$total <- intercepted.seed.wide$Not_intercepted + intercepted.seed.wide$intercepted
intercepted.seed.wide$prop.intercepted <- intercepted.seed.wide$intercepted/intercepted.seed.wide$total
prop.intercepted.seed <- paste0(intercepted.seed.wide$intercepted, "/", intercepted.seed.wide$total)

p_interceptedseed <- ggplot(data=intercepted.seed.wide, aes(x=Seed_Dispersal, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "***")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.intercepted.seed), color="white")

p_intercepted <- plot_grid(p_interceptedefn, p_intercepteddom, p_interceptedseed, nrow=1, labels="AUTO")

binomial9 <- glmer(InterceptedY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|genus), data=subset(area, introducedY == 1), family="binomial")
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## unable to evaluate scaled gradient

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge: degenerate Hessian with 1 negative eigenvalues

``` r
summary(binomial9)
```

    ## Warning in vcov.merMod(object, use.hessian = use.hessian): variance-covariance matrix computed from finite-difference Hessian is
    ## not positive definite or contains NA values: falling back to var-cov estimated from RX

    ## Warning in vcov.merMod(object, correlation = correlation, sigm = sig): variance-covariance matrix computed from finite-difference Hessian is
    ## not positive definite or contains NA values: falling back to var-cov estimated from RX

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## InterceptedY ~ EFN + Seed_Dispersal + Domatia + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | genus)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    435.6    462.4   -210.8    421.6      334 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9479 -0.8708  0.4202  0.7912  1.6650 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 0.3949   0.6284  
    ## Number of obs: 341, groups:  genus, 69
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                -0.1730     0.1743  -0.992 0.321069    
    ## EFN1                        0.5089     0.4411   1.154 0.248608    
    ## Seed_Dispersal1             1.2026     0.3292   3.653 0.000259 ***
    ## Domatia1                   15.5456  1459.5820   0.011 0.991502    
    ## scale(abs_lat_native)      -0.3473     0.1291  -2.690 0.007154 ** 
    ## scale(total.area.native)    0.1422     0.1317   1.079 0.280566    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Sd_Ds1 Domat1 sc(__)
    ## EFN1        -0.151                            
    ## Sed_Dsprsl1 -0.379 -0.148                     
    ## Domatia1     0.000  0.000  0.000              
    ## scl(bs_lt_)  0.054  0.085 -0.267  0.000       
    ## scl(ttl.r.)  0.067 -0.122 -0.259  0.000 -0.056
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## unable to evaluate scaled gradient
    ## Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

``` r
Anova(binomial9, type=3)
```

    ## Warning in vcov.merMod(mod, complete = FALSE): variance-covariance matrix computed from finite-difference Hessian is
    ## not positive definite or contains NA values: falling back to var-cov estimated from RX

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: InterceptedY
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)               0.9846  1  0.3210692    
    ## EFN                       1.3311  1  0.2486081    
    ## Seed_Dispersal           13.3431  1  0.0002594 ***
    ## Domatia                   0.0001  1  0.9915021    
    ## scale(abs_lat_native)     7.2339  1  0.0071540 ** 
    ## scale(total.area.native)  1.1643  1  0.2805660    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial9)
```

![](README_files/figure-gfm/Introduction%20mode-3.png)<!-- -->

``` r
fig7 <- plot_grid(p_intercepted, p_indoor, p_exotic, nrow=3)
fig7
```

![](README_files/figure-gfm/Introduction%20mode-4.png)<!-- -->

``` r
save_plot("Figure7.pdf", fig7, base_height=8, base_width=8)
```

### PGLS Models

``` r
#all ant mutualisms
#Native range size
ant_tree <- read.tree("Ant_tree.tre") #Reading in ant phylogeny
phy_int <- intersect(ant_tree$tip.label, area$Phy)
phy_diff <- setdiff(ant_tree$tip.label, area$Phy)
pruned_ant_tree <- drop.tip(ant_tree, as.character(phy_diff)) #pruning tree to contain only tips in the dataset

area_phy <- area[area$Phy %in% phy_int, c("Phy", "abs_lat_native", "total.area.introduced", "total.area.native", "EFN", "Domatia", "Seed_Dispersal")] #dataset for pgls
area_phy <- area_phy[complete.cases(area_phy), ]
#write.csv(efn_area_phy, "Ant_list_EFN.csv")
#efn_area_phy <- efn_area_phy[match(pruned_ant_tree$tip.label, efn_area_phy$Phy),]

#PGLS model
#a1 <- gls(log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + EFN*Domatia + Domatia*Seed_Dispersal + Seed_Dispersal*EFN + abs_lat_native, 
                #correlation = corPagel(1, phy = pruned_ant_tree, form = ~ Phy), 
                #method = "ML", data = area_phy) 
#summary(a1) 

#removing interactions
a11 <- gls(log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + abs_lat_native, 
                correlation = corPagel(1, phy = pruned_ant_tree, form = ~ Phy), 
                method = "ML", data = area_phy) 
summary(a11)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + abs_lat_native 
    ##   Data: area_phy 
    ##       AIC      BIC    logLik
    ##   2607.99 2640.506 -1296.995
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.5359385 
    ## 
    ## Coefficients:
    ##                     Value Std.Error  t-value p-value
    ## (Intercept)     14.366824 0.4167675 34.47204  0.0000
    ## EFN1             0.614639 0.1743272  3.52578  0.0004
    ## Domatia1        -0.042505 0.2586044 -0.16436  0.8695
    ## Seed_Dispersal1  0.932888 0.1309378  7.12466  0.0000
    ## abs_lat_native   0.004148 0.0043174  0.96071  0.3370
    ## 
    ##  Correlation: 
    ##                 (Intr) EFN1   Domat1 Sd_Ds1
    ## EFN1            -0.026                     
    ## Domatia1        -0.029 -0.138              
    ## Seed_Dispersal1 -0.015 -0.185  0.032       
    ## abs_lat_native  -0.163  0.028  0.068 -0.026
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -5.0299546 -0.1284837  0.4385483  0.8795089  1.9368672 
    ## 
    ## Residual standard error: 1.591093 
    ## Degrees of freedom: 769 total; 764 residual

``` r
#Diagnostic plots
plot(a11$residuals, a11$fitted)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-1.png)<!-- -->

``` r
qqnorm(a11$residuals)
qqline(a11$residuals)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-2.png)<!-- -->

``` r
anova(a11)
```

    ## Denom. DF: 764 
    ##                numDF   F-value p-value
    ## (Intercept)        1 1252.3710  <.0001
    ## EFN                1   23.8944  <.0001
    ## Domatia            1    0.2217  0.6379
    ## Seed_Dispersal     1   51.1537  <.0001
    ## abs_lat_native     1    0.9230  0.3370

``` r
#Introduced range size 
a2 <- gls(log(total.area.introduced+1) ~  EFN + Domatia + Seed_Dispersal +  total.area.native + abs_lat_native,  correlation = corPagel(1, phy = pruned_ant_tree, form = ~Phy), method = "ML", data = area_phy)  
summary(a2)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(total.area.introduced + 1) ~ EFN + Domatia + Seed_Dispersal +      total.area.native + abs_lat_native 
    ##   Data: area_phy 
    ##        AIC      BIC    logLik
    ##   4580.651 4617.812 -2282.326
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.5040291 
    ## 
    ## Coefficients:
    ##                        Value Std.Error   t-value p-value
    ## (Intercept)        1.6508696 1.4447771  1.142647  0.2535
    ## EFN1               1.9445287 0.6424026  3.026963  0.0026
    ## Domatia1          -0.6669233 0.9330685 -0.714764  0.4750
    ## Seed_Dispersal1    2.4630599 0.4944578  4.981335  0.0000
    ## total.area.native  0.0000002 0.0000000  7.260721  0.0000
    ## abs_lat_native    -0.0693057 0.0155102 -4.468397  0.0000
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN1   Domat1 Sd_Ds1 ttl.r.
    ## EFN1              -0.013                            
    ## Domatia1          -0.032 -0.141                     
    ## Seed_Dispersal1    0.005 -0.121  0.025              
    ## total.area.native -0.069 -0.192  0.022 -0.288       
    ## abs_lat_native    -0.166  0.038  0.070 -0.014 -0.048
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -1.6588226 -0.3843118 -0.1282512  0.1180702  2.8455384 
    ## 
    ## Residual standard error: 5.628969 
    ## Degrees of freedom: 769 total; 763 residual

``` r
## Diagnostic plots
plot(a2$residuals, a2$fitted)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-3.png)<!-- -->

``` r
qqnorm(a2$residuals)
qqline(a2$residuals)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-4.png)<!-- -->
