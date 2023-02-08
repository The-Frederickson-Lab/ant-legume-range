Generalized mutualism promotes range expansion in both ant and plant
partners
================
Pooja Nathan and Megan Frederickson
2023-02-07

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
library(corrplot)
library(Hmisc)
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
range <- read.csv("inv_dat_by_species_simonsen2017.csv") #Read in legume range data
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

mycorrhizae.sum <- mycorrhizae %>% group_by(Phy, species, RangeY, domY, efnY) %>% dplyr::summarize(n=n(), sum.AM = sum(AM), sum.EM=sum(EM))
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
mycorrhizae <- mycorrhizae %>% group_by(Phy, RangeY, domY, efnY) %>% dplyr::summarize(n.records=sum(n), sum.AM = sum(sum.AM), sum.EM=sum(sum.EM))
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
legume_range_df$woody <- as.numeric(legume_range_df$woody)

##Collapse all mycorrhizal fungi types into a single yes/no category
legume_range_df$myco <- ifelse(legume_range_df$AM == "Y" | legume_range_df$EM == "Y", 1, ifelse(legume_range_df$AM == "N" & legume_range_df$EM == "N", 0, NA))
legume_range_df$myco <- as.factor(legume_range_df$myco)

df <- legume_range_df
#annual and woody are highly correlated
#combining them into a single variable called Natural History
#df$nathis <- df$annual + df$woody
#df$nathis <- as.factor(df$nathis)

summary.efn <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.efn)
```

| EFN |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|-----:|--------------------:|------------------:|------------------:|
| 0   | 3697 |            1.365702 |          5.330667 |         0.0876712 |
| 1   |  280 |            5.257143 |         10.864345 |         0.6492688 |

How many legumes with vs. without domatia do we have range size data
for?

``` r
summary.dom <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.dom)
```

| Domatia |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:--------|-----:|--------------------:|------------------:|------------------:|
| 0       | 3953 |           1.6463445 |          5.990833 |         0.0952848 |
| 1       |   24 |           0.5416667 |          1.178767 |         0.2406149 |

How many legumes that do vs. do not form nodules do we have range size
data for?

``` r
summary.fix <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.fix)
```

| fixer |    n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:------|-----:|--------------------:|------------------:|------------------:|
| 0     |  396 |            2.482323 |          6.839364 |         0.3436910 |
| 1     | 3581 |            1.546495 |          5.864139 |         0.0979946 |

How many legumes do vs. do not associate with mycorrhizae do we have
range size data for?

``` r
summary.myco <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(myco) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.myco)
```

| myco |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:-----|----:|--------------------:|------------------:|------------------:|
| 0    |  33 |            3.636364 |          7.192957 |         1.2521332 |
| 1    | 690 |            5.998551 |         11.970876 |         0.4557235 |

``` r
summary.myco2 <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(AM,EM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
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
summary.efn <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.efn)
```

| EFN |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| 0   | 716 |            7.051676 |          10.33137 |         0.3861014 |
| 1   | 135 |           10.903704 |          13.55468 |         1.1666015 |

``` r
p_EFN <- ggplot(data=summary.efn, aes(x=EFN, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "***")

summary.efn2 <- ungroup(df %>% group_by(EFN, introducedY) %>% dplyr::summarize(n=n()))
summary.efn2.wide <- spread(summary.efn2, key = introducedY, value=n)
colnames(summary.efn2.wide) <- c("EFN","Not_introduced",  "Introduced")
summary.efn2.wide$total <- summary.efn2.wide$Not_introduced+summary.efn2.wide$Introduced
summary.efn2.wide$prop.introduced <- summary.efn2.wide$Introduced/(summary.efn2.wide$total)
prop.efn <- paste0(summary.efn2.wide$Introduced, "/", summary.efn2.wide$total)

inset_p_EFN <- ggplot(data=summary.efn2.wide, aes(x=EFN, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=EFN, y=0.05, label=prop.efn), color="white")

#Domatia figure
summary.dom <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.dom)
```

| Domatia |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:--------|----:|--------------------:|------------------:|------------------:|
| 0       | 844 |            7.710901 |         11.019742 |         0.3793152 |
| 1       |   7 |            1.857143 |          1.573592 |         0.5947617 |

``` r
p_dom <- ggplot(data=summary.dom, aes(x=Domatia, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "**")

summary.dom2 <- ungroup(df %>% group_by(Domatia, introducedY) %>% dplyr::summarize(n=n()))
summary.dom2.wide <- spread(summary.dom2, key = introducedY, value=n)
colnames(summary.dom2.wide) <- c("Domatia","Not_introduced",  "Introduced")
summary.dom2.wide$total <- summary.dom2.wide$Not_introduced+summary.dom2.wide$Introduced
summary.dom2.wide$prop.introduced <- summary.dom2.wide$Introduced/summary.dom2.wide$total
prop.dom <- paste0(summary.dom2.wide$Introduced, "/", summary.dom2.wide$total)

inset_p_dom <- ggplot(data=summary.dom2.wide, aes(x=Domatia, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.dom), color="white")

#Nodules figure
summary.fix <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.fix)
```

| fixer |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:------|----:|--------------------:|------------------:|------------------:|
| 0     | 103 |            9.543689 |          10.63455 |         1.0478534 |
| 1     | 748 |            7.403743 |          11.01733 |         0.4028336 |

``` r
p_fix <- ggplot(data=summary.fix, aes(x=fixer, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "ns")

summary.fix2 <- ungroup(df %>% group_by(fixer, introducedY) %>% dplyr::summarize(n=n()))
summary.fix2.wide <- spread(summary.fix2, key = introducedY, value=n)
colnames(summary.fix2.wide) <- c("fixer","Not_introduced",  "Introduced")
summary.fix2.wide$total <- summary.fix2.wide$Not_introduced+summary.fix2.wide$Introduced
summary.fix2.wide$prop.introduced <- summary.fix2.wide$Introduced/summary.fix2.wide$total
prop.fix <- paste0(summary.fix2.wide$Introduced, "/", summary.fix2.wide$total)
  
inset_p_fix <- ggplot(data=summary.fix2.wide, aes(x=fixer, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Nodules")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=fixer, y=0.05, label=prop.fix), color="white")

#Mycorrhizae figure
summary.AM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(AM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.AM)
```

| AM  |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| N   |  22 |            6.954546 |          8.126836 |         1.7326471 |
| Y   | 308 |           13.331169 |         14.933754 |         0.8509296 |

``` r
p_AM <- ggplot(data=summary.AM, aes(x=AM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns")

summary.EM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(EM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.EM)
```

| EM  |   n | mean_num_introduced | sd_num_introduced | se_num_introduced |
|:----|----:|--------------------:|------------------:|------------------:|
| N   | 295 |           13.264407 |          14.92901 |         0.8692007 |
| Y   |  35 |            9.885714 |          11.88863 |         2.0095451 |

``` r
p_EM <- ggplot(data=summary.EM, aes(x=EM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("EM")+geom_text(aes(x=EM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns")

#Or
summary.AMEM2 <- ungroup(subset(df, !is.na(myco)) %>% group_by(AM, EM, introducedY) %>% dplyr::summarize(n=n()))
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

summary.efn.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))

p_EFN_area <- ggplot(data=summary.efn.area, aes(x=EFN, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.dom.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.dom.area)
```

| Domatia |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:--------|-----:|---------------------:|-------------------:|-------------------:|
| 0       | 3953 |         1.582357e+12 |       6.502144e+12 |       103417297082 |
| 1       |   24 |         5.290977e+11 |       1.104627e+12 |       225481104596 |

``` r
p_dom_area <- ggplot(data=summary.dom.area, aes(x=Domatia, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.fix.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.fix.area)
```

| fixer |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:------|-----:|---------------------:|-------------------:|-------------------:|
| 0     |  396 |         1.995301e+12 |       6.040271e+12 |       303535045685 |
| 1     | 3581 |         1.529633e+12 |       6.529860e+12 |       109119327208 |

``` r
p_fix_area <- ggplot(data=summary.fix.area, aes(x=fixer, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range area (sq. km)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.myco.area <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(myco) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.myco.area)
```

| myco |    n | mean_area_introduced | sd_area_introduced | se_area_introduced |
|:-----|-----:|---------------------:|-------------------:|-------------------:|
| 0    |   33 |         3.194037e+12 |       8.008639e+12 |       1.394125e+12 |
| 1    |  690 |         5.765715e+12 |       1.325680e+13 |       5.046778e+11 |
| NA   | 3254 |         6.711765e+11 |       2.992716e+12 |       5.246343e+10 |

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

summary.efn.dom.native <- ungroup(df %>% group_by(EFN, Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_EFN_dom_native <- ggplot(data=summary.efn.dom.native, aes(x=EFN, y=mean_native, color=Domatia))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native, color=Domatia), width=er_width)+ geom_line(aes(group=Domatia), linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("EFNs")+geom_text(aes(x=c(0.8,1.2,1.8,2.2), y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+scale_color_grey(labels=c("No", "Yes"))+theme(legend.position = c(0.1, 0.8))+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")+annotate("text", x = 1.5, y = 3.9e+12, label = "***")

summary.dom.native <- ungroup(df %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_dom_native <- ggplot(data=summary.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)

summary.fix.native <- ungroup(df %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_fix_native <- ggplot(data=summary.fix.native, aes(x=fixer, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

summary.AM.native <- ungroup(df %>% group_by(AM) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.AM.native)
```

| AM  |    n |  mean_native |    sd_native |    se_native |
|:----|-----:|-------------:|-------------:|-------------:|
| N   |   71 | 6.641775e+12 | 7.620421e+12 | 904377608987 |
| Y   |  652 | 9.645563e+12 | 1.042585e+13 | 408307721818 |
| NA  | 3254 | 4.185009e+12 | 4.948411e+12 |  86747509406 |

``` r
p_AM_native <- ggplot(data=subset(summary.AM.native, !is.na(AM)), aes(x=AM, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

summary.EM.native <- ungroup(df %>% group_by(EM) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.EM.native)
```

| EM  |    n |  mean_native |    sd_native |    se_native |
|:----|-----:|-------------:|-------------:|-------------:|
| N   |  605 | 1.022114e+13 | 1.082192e+13 | 439973745086 |
| Y   |  118 | 4.887136e+12 | 3.986129e+12 | 366952887362 |
| NA  | 3254 | 4.185009e+12 | 4.948411e+12 |  86747509406 |

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
colnames(range_pgls)
```

    ##  [1] "Phy"                   "introduced"            "native"               
    ##  [4] "num_introduced"        "total_area_introduced" "total_area_native"    
    ##  [7] "lat_native"            "tribe"                 "genus"                
    ## [10] "subfamily"             "tribe_ncbi"            "fixer"                
    ## [13] "abs_lat_native"        "woody"                 "annual"               
    ## [16] "uses_num_uses"         "submitted_name"        "matched_name"         
    ## [19] "data_source_title"     "score"                 "num_words"            
    ## [22] "diff"                  "RangeY"                "domY"                 
    ## [25] "efnY"                  "n.records"             "sum.AM"               
    ## [28] "sum.EM"                "AM"                    "EM"                   
    ## [31] "EFN"                   "matched_name_EFN"      "Domatia"              
    ## [34] "myco"                  "introducedY"           "Phy2"

``` r
which(colSums(is.na(range_pgls))>0) #Check which columns have NAs
```

    ##       lat_native   abs_lat_native           RangeY             domY 
    ##                7               13               23               24 
    ##             efnY        n.records           sum.AM           sum.EM 
    ##               25               26               27               28 
    ##               AM               EM matched_name_EFN             myco 
    ##               29               30               32               34

``` r
range_pgls <-range_pgls[,-c(23:30, 32, 34)] #Remove some unneeded columns
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

pgls1 <- gls(log(num_introduced + 1) ~  EFN+fixer+Domatia+total_area_native+abs_lat_native+uses_num_uses+woody+annual, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method="ML", data = range_pgls) 
summary(pgls1)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(num_introduced + 1) ~ EFN + fixer + Domatia + total_area_native +      abs_lat_native + uses_num_uses + woody + annual 
    ##   Data: range_pgls 
    ##        AIC      BIC    logLik
    ##   2679.213 2735.385 -1328.606
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.2976025 
    ## 
    ## Coefficients:
    ##                        Value  Std.Error  t-value p-value
    ## (Intercept)        0.1447488 0.19825832  0.73010  0.4655
    ## EFN1               0.4658576 0.07389168  6.30460  0.0000
    ## fixer1             0.1076161 0.11480643  0.93737  0.3488
    ## Domatia1          -0.2740403 0.31049160 -0.88260  0.3776
    ## total_area_native  0.0000000 0.00000000 -1.76214  0.0783
    ## abs_lat_native    -0.0024371 0.00204928 -1.18925  0.2346
    ## uses_num_uses      0.4164800 0.01247662 33.38083  0.0000
    ## woody              0.0908155 0.06929044  1.31065  0.1902
    ## annual             0.0902686 0.07223048  1.24973  0.2116
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN1   fixer1 Domat1 ttl_r_ abs_l_ uss_n_ woody 
    ## EFN1               0.006                                                 
    ## fixer1            -0.232 -0.069                                          
    ## Domatia1           0.003 -0.001 -0.097                                   
    ## total_area_native -0.097 -0.021 -0.038 -0.006                            
    ## abs_lat_native    -0.239  0.024 -0.054  0.053  0.090                     
    ## uses_num_uses     -0.049 -0.118  0.068  0.020 -0.429  0.020              
    ## woody             -0.348 -0.004  0.021 -0.004  0.109  0.115 -0.095       
    ## annual            -0.165 -0.009 -0.014  0.006  0.113  0.126  0.008  0.380
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -3.5168049 -0.3719842 -0.2767498  0.3861284  3.2323950 
    ## 
    ## Residual standard error: 0.8140867 
    ## Degrees of freedom: 1220 total; 1211 residual

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
            total_area_native + abs_lat_native + uses_num_uses+ annual + woody, 
            correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), 
            method = "ML", data = range_pgls) 
summary(pgls2)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(total_area_introduced/1e+06 + 1) ~ EFN + Domatia + fixer +      total_area_native + abs_lat_native + uses_num_uses + annual +      woody 
    ##   Data: range_pgls 
    ##        AIC      BIC    logLik
    ##   7636.922 7693.094 -3807.461
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.1899116 
    ## 
    ## Coefficients:
    ##                       Value Std.Error   t-value p-value
    ## (Intercept)        0.452192 1.2753895  0.354552  0.7230
    ## EFN1               3.171629 0.5659303  5.604275  0.0000
    ## Domatia1          -1.320612 2.3389460 -0.564619  0.5724
    ## fixer1             0.715989 0.8145744  0.878973  0.3796
    ## total_area_native  0.000000 0.0000000  0.522089  0.6017
    ## abs_lat_native     0.026247 0.0151932  1.727547  0.0843
    ## uses_num_uses      2.232167 0.0951167 23.467655  0.0000
    ## annual             1.370377 0.5497385  2.492780  0.0128
    ## woody              1.323314 0.5166018  2.561574  0.0105
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN1   Domat1 fixer1 ttl_r_ abs_l_ uss_n_ annual
    ## EFN1               0.007                                                 
    ## Domatia1           0.002  0.003                                          
    ## fixer1            -0.282 -0.071 -0.096                                   
    ## total_area_native -0.118 -0.025 -0.005 -0.041                            
    ## abs_lat_native    -0.277  0.025  0.060 -0.068  0.097                     
    ## uses_num_uses     -0.060 -0.116  0.021  0.078 -0.430  0.018              
    ## annual            -0.197 -0.010  0.007 -0.017  0.111  0.129  0.009       
    ## woody             -0.405 -0.003 -0.005  0.027  0.120  0.130 -0.093  0.395
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.8746120 -0.5982965 -0.4455098  0.7294161  2.5320401 
    ## 
    ## Residual standard error: 5.884291 
    ## Degrees of freedom: 1220 total; 1211 residual

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
pgls3 <- gls(log((total_area_native/1e+6) + 1) ~ EFN*Domatia + fixer+ abs_lat_native+ annual + woody + uses_num_uses, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method = "ML", data = range_pgls)
summary(pgls3)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log((total_area_native/1e+06) + 1) ~ EFN * Domatia + fixer +      abs_lat_native + annual + woody + uses_num_uses 
    ##   Data: range_pgls 
    ##        AIC      BIC    logLik
    ##   4095.409 4151.582 -2036.704
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy2 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.4607838 
    ## 
    ## Coefficients:
    ##                    Value Std.Error  t-value p-value
    ## (Intercept)    15.269905 0.4439655 34.39435  0.0000
    ## EFN1           -0.078439 0.1309877 -0.59883  0.5494
    ## Domatia1        0.515727 0.6319693  0.81606  0.4146
    ## fixer1          0.040534 0.2212697  0.18319  0.8547
    ## abs_lat_native -0.008597 0.0037595 -2.28685  0.0224
    ## annual         -0.178108 0.1280049 -1.39142  0.1644
    ## woody          -0.409857 0.1267516 -3.23355  0.0013
    ## uses_num_uses   0.250756 0.0201408 12.45014  0.0000
    ## EFN1:Domatia1  -3.307219 1.3610321 -2.42993  0.0152
    ## 
    ##  Correlation: 
    ##                (Intr) EFN1   Domat1 fixer1 abs_l_ annual woody  uss_n_
    ## EFN1            0.005                                                 
    ## Domatia1        0.004  0.032                                          
    ## fixer1         -0.185 -0.070 -0.096                                   
    ## abs_lat_native -0.190  0.026  0.041 -0.042                            
    ## annual         -0.122 -0.005  0.005 -0.008  0.113                     
    ## woody          -0.278 -0.005 -0.002  0.020  0.094  0.355              
    ## uses_num_uses  -0.076 -0.145  0.000  0.047  0.064  0.064 -0.062       
    ## EFN1:Domatia1  -0.004 -0.080 -0.456  0.039 -0.005  0.000 -0.001  0.029
    ## 
    ## Standardized residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -7.48404947 -0.37504275  0.09830079  0.54896201  1.53698574 
    ## 
    ## Residual standard error: 1.609154 
    ## Degrees of freedom: 1220 total; 1211 residual

``` r
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

    ##       lat_native   abs_lat_native matched_name_EFN 
    ##                7               13               32

``` r
range_myco_pgls <-range_myco_pgls[,-c(32)] #Remove some unneeded columns
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
#Checking correlation between variables
df_cor <- data.frame(as.numeric(df$fixer), as.numeric(df$EFN), as.numeric(df$Domatia), as.numeric(df$EM), as.numeric(df$AM), #mutualisms
                     df$woody, df$annual, df$uses_num_uses, df$abs_lat_native) #covariates
# <- cor(df_cor)
colnames(df_cor) <- c("fixer", "EFN", "Domatia", "EM", "AM", "Woody", "Annual", "Uses", "Abs_lat")
corr <- rcorr(as.matrix(df_cor))

plotr <- corrplot(corr$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) #correlation coeffs
```

![](README_files/figure-gfm/Legume%20mixed%20models-1.png)<!-- -->

``` r
plotp <- corrplot(corr$P, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) #p-values
```

![](README_files/figure-gfm/Legume%20mixed%20models-2.png)<!-- -->

``` r
#Fit binomial model for whether or not a legume species has been introduced
binomial1 <- glmer(introducedY~fixer*uses_num_uses + fixer + EFN+Domatia+scale(total_area_native)+woody+annual+scale(abs_lat_native)+(1|tribe_ncbi), data=df, family="binomial",glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(binomial1) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: introducedY ~ fixer * uses_num_uses + fixer + EFN + Domatia +  
    ##     scale(total_area_native) + woody + annual + scale(abs_lat_native) +  
    ##     (1 | tribe_ncbi)
    ##    Data: df
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2544.9   2612.8  -1261.4   2522.9     3523 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -14.7388  -0.3819  -0.2651  -0.1467   6.2195 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.5034   0.7095  
    ## Number of obs: 3534, groups:  tribe_ncbi, 51
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -3.20650    0.31480 -10.186  < 2e-16 ***
    ## fixer1                    0.02091    0.27184   0.077 0.938700    
    ## uses_num_uses             0.75967    0.10396   7.307 2.73e-13 ***
    ## EFN1                      1.25903    0.19517   6.451 1.11e-10 ***
    ## Domatia1                  0.39750    0.64317   0.618 0.536553    
    ## scale(total_area_native)  0.17335    0.05770   3.004 0.002663 ** 
    ## woody                     0.78305    0.16288   4.808 1.53e-06 ***
    ## annual                    0.68973    0.18810   3.667 0.000246 ***
    ## scale(abs_lat_native)     0.09114    0.06933   1.315 0.188654    
    ## fixer1:uses_num_uses      0.29484    0.11605   2.541 0.011063 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) fixer1 uss_n_ EFN1   Domat1 scl(t__) woody  annual scl(b__)
    ## fixer1      -0.761                                                            
    ## uses_num_ss -0.464  0.521                                                     
    ## EFN1        -0.060  0.041  0.015                                              
    ## Domatia1     0.027 -0.063  0.033 -0.013                                       
    ## scl(ttl_r_) -0.046  0.032 -0.055 -0.033  0.006                                
    ## woody       -0.462  0.056 -0.015 -0.034 -0.001  0.142                         
    ## annual      -0.195 -0.057 -0.009 -0.018  0.020  0.079    0.419                
    ## scl(bs_lt_)  0.018 -0.104  0.003  0.027  0.095  0.129    0.150  0.130         
    ## fxr1:ss_nm_  0.403 -0.543 -0.884 -0.001 -0.017 -0.109    0.015  0.056 -0.012

``` r
Anova(binomial1, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)              103.7524  1  < 2.2e-16 ***
    ## fixer                      0.0059  1  0.9387000    
    ## uses_num_uses             53.3929  1  2.731e-13 ***
    ## EFN                       41.6145  1  1.112e-10 ***
    ## Domatia                    0.3820  1  0.5365527    
    ## scale(total_area_native)   9.0249  1  0.0026632 ** 
    ## woody                     23.1122  1  1.528e-06 ***
    ## annual                    13.4461  1  0.0002455 ***
    ## scale(abs_lat_native)      1.7281  1  0.1886537    
    ## fixer:uses_num_uses        6.4551  1  0.0110633 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial1)
```

![](README_files/figure-gfm/Legume%20mixed%20models-3.png)<!-- -->

``` r
#EFN, Domatia interaction with human uses was non-significant

#Fit binomial model without random effect
binomial2 <- glm(introducedY~fixer*uses_num_uses+EFN*Domatia+fixer+scale(total_area_native)+woody+annual+scale(abs_lat_native)+uses_num_uses, data=df, family="binomial")
summary(binomial2) 
```

    ## 
    ## Call:
    ## glm(formula = introducedY ~ fixer * uses_num_uses + EFN * Domatia + 
    ##     fixer + scale(total_area_native) + woody + annual + scale(abs_lat_native) + 
    ##     uses_num_uses, family = "binomial", data = df)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.3161  -0.5035  -0.4233  -0.3193   2.5368  
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.58318    0.24473 -10.555  < 2e-16 ***
    ## fixer1                   -0.22870    0.22398  -1.021 0.307206    
    ## uses_num_uses             0.73919    0.09651   7.659 1.87e-14 ***
    ## EFN1                      1.28905    0.17224   7.484 7.21e-14 ***
    ## Domatia1                  0.46998    0.75473   0.623 0.533474    
    ## scale(total_area_native)  0.20525    0.05473   3.750 0.000176 ***
    ## woody                     0.52827    0.14376   3.675 0.000238 ***
    ## annual                    0.98970    0.16762   5.904 3.54e-09 ***
    ## scale(abs_lat_native)     0.13697    0.05312   2.579 0.009918 ** 
    ## fixer1:uses_num_uses      0.30376    0.10858   2.798 0.005149 ** 
    ## EFN1:Domatia1             0.35515    1.58532   0.224 0.822740    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3814.5  on 3533  degrees of freedom
    ## Residual deviance: 2636.7  on 3523  degrees of freedom
    ##   (443 observations deleted due to missingness)
    ## AIC: 2658.7
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
Anova(binomial2, type=3)
```

    ## Analysis of Deviance Table (Type III tests)
    ## 
    ## Response: introducedY
    ##                          LR Chisq Df Pr(>Chisq)    
    ## fixer                       1.007  1  0.3157180    
    ## uses_num_uses              90.863  1  < 2.2e-16 ***
    ## EFN                        53.176  1  3.050e-13 ***
    ## Domatia                     0.348  1  0.5554063    
    ## scale(total_area_native)   14.225  1  0.0001622 ***
    ## woody                      14.038  1  0.0001792 ***
    ## annual                     34.796  1  3.662e-09 ***
    ## scale(abs_lat_native)       6.688  1  0.0097076 ** 
    ## fixer:uses_num_uses         6.999  1  0.0081575 ** 
    ## EFN:Domatia                 0.051  1  0.8215439    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-4.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-5.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-6.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-7.png)<!-- -->

``` r
#Compare models
anova(binomial1, binomial2) #Random effect improves model fit to data
```

    ## Data: df
    ## Models:
    ## binomial1: introducedY ~ fixer * uses_num_uses + fixer + EFN + Domatia + scale(total_area_native) + woody + annual + scale(abs_lat_native) + (1 | tribe_ncbi)
    ## binomial2: introducedY ~ fixer * uses_num_uses + EFN * Domatia + fixer + scale(total_area_native) + woody + annual + scale(abs_lat_native) + uses_num_uses
    ##           npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
    ## binomial1   11 2544.9 2612.8 -1261.4   2522.9                    
    ## binomial2   11 2658.7 2726.6 -1318.4   2636.7     0  0

``` r
#Fit a linear mixed model for how many new ranges an introduced legume has established in, and how much total area they cover
legume_range_df_introducedY <- subset(df, num_introduced >0) #Filter to species with 1+ introduced ranges

#Number of introduced ranges
lmer1 <- lmer(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+woody+uses_num_uses+(1|tribe_ncbi), data=legume_range_df_introducedY)
summary(lmer1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) +  
    ##     scale(total_area_native) + annual + woody + uses_num_uses +  
    ##     (1 | tribe_ncbi)
    ##    Data: legume_range_df_introducedY
    ## 
    ## REML criterion at convergence: 2208.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.09452 -0.76189 -0.02237  0.68385  2.70244 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.04084  0.2021  
    ##  Residual               0.83733  0.9151  
    ## Number of obs: 814, groups:  tribe_ncbi, 39
    ## 
    ## Fixed effects:
    ##                           Estimate Std. Error t value
    ## (Intercept)               0.852767   0.162638   5.243
    ## EFN1                      0.343026   0.098576   3.480
    ## Domatia1                 -0.986516   0.383716  -2.571
    ## fixer1                   -0.210014   0.136135  -1.543
    ## scale(abs_lat_native)    -0.064120   0.042109  -1.523
    ## scale(total_area_native) -0.051380   0.037046  -1.387
    ## annual                    0.007123   0.119607   0.060
    ## woody                    -0.223967   0.099458  -2.252
    ## uses_num_uses             0.321275   0.016440  19.542
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 scl(b__) scl(t__) annual woody 
    ## EFN1        -0.049                                                     
    ## Domatia1     0.025 -0.013                                              
    ## fixer1      -0.780  0.031 -0.040                                       
    ## scl(bs_lt_)  0.022  0.062  0.064 -0.131                                
    ## scl(ttl_r_)  0.078 -0.029 -0.004 -0.109  0.061                         
    ## annual      -0.256  0.013 -0.002 -0.029 -0.008    0.187                
    ## woody       -0.543 -0.085 -0.012  0.118  0.135    0.168    0.418       
    ## uses_num_ss -0.284 -0.044 -0.019  0.087  0.036   -0.376    0.012 -0.053

``` r
Anova(lmer1, type=3,)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               27.4926  1  1.577e-07 ***
    ## EFN                       12.1091  1  0.0005018 ***
    ## Domatia                    6.6098  1  0.0101419 *  
    ## fixer                      2.3799  1  0.1229066    
    ## scale(abs_lat_native)      2.3187  1  0.1278244    
    ## scale(total_area_native)   1.9236  1  0.1654623    
    ## annual                     0.0035  1  0.9525140    
    ## woody                      5.0710  1  0.0243295 *  
    ## uses_num_uses            381.8802  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer1)
```

![](README_files/figure-gfm/Legume%20mixed%20models-8.png)<!-- -->

``` r
#Removed the interaction when nonsig or the number of species showing trait combinations was too small

#Number of introduced ranges without random effect of genus
lm2 <- lm(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+woody+uses_num_uses, data=legume_range_df_introducedY)
summary(lm2)
```

    ## 
    ## Call:
    ## lm(formula = log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + 
    ##     scale(total_area_native) + annual + woody + uses_num_uses, 
    ##     data = legume_range_df_introducedY)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0074 -0.7321 -0.0230  0.6482  2.5258 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)               0.96476    0.13783   7.000 5.41e-12 ***
    ## EFN1                      0.36826    0.09165   4.018 6.42e-05 ***
    ## Domatia1                 -1.03072    0.38320  -2.690  0.00730 ** 
    ## fixer1                   -0.27019    0.10700  -2.525  0.01176 *  
    ## scale(abs_lat_native)    -0.04390    0.03554  -1.235  0.21702    
    ## scale(total_area_native) -0.04904    0.03648  -1.344  0.17921    
    ## annual                    0.07310    0.11491   0.636  0.52487    
    ## woody                    -0.28276    0.09328  -3.031  0.00251 ** 
    ## uses_num_uses             0.32206    0.01644  19.587  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.9297 on 805 degrees of freedom
    ##   (37 observations deleted due to missingness)
    ## Multiple R-squared:  0.3763, Adjusted R-squared:  0.3701 
    ## F-statistic: 60.71 on 8 and 805 DF,  p-value: < 2.2e-16

``` r
Anova(lm2, type=3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(num_introduced)
    ##                          Sum Sq  Df  F value    Pr(>F)    
    ## (Intercept)               42.35   1  48.9949  5.41e-12 ***
    ## EFN                       13.95   1  16.1441  6.42e-05 ***
    ## Domatia                    6.25   1   7.2349  0.007298 ** 
    ## fixer                      5.51   1   6.3757  0.011760 *  
    ## scale(abs_lat_native)      1.32   1   1.5264  0.217019    
    ## scale(total_area_native)   1.56   1   1.8073  0.179214    
    ## annual                     0.35   1   0.4047  0.524868    
    ## woody                      7.94   1   9.1885  0.002513 ** 
    ## uses_num_uses            331.62   1 383.6613 < 2.2e-16 ***
    ## Residuals                695.81 805                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lm2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-9.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-10.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-11.png)<!-- -->![](README_files/figure-gfm/Legume%20mixed%20models-12.png)<!-- -->

``` r
#Compare models
anova(lmer1, lm2) #Random effect improves model fit to data
```

    ## Data: legume_range_df_introducedY
    ## Models:
    ## lm2: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + scale(total_area_native) + annual + woody + uses_num_uses
    ## lmer1: log(num_introduced) ~ EFN + Domatia + fixer + scale(abs_lat_native) + scale(total_area_native) + annual + woody + uses_num_uses + (1 | tribe_ncbi)
    ##       npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
    ## lm2     10 2202.3 2249.3 -1091.2   2182.3                       
    ## lmer1   11 2199.8 2251.5 -1088.9   2177.8 4.5284  1    0.03334 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#Introduced area (if introduced at all)
lmer2 <- lmer(log(total_area_introduced/1e+6)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+woody+uses_num_uses+(1|tribe_ncbi), data=legume_range_df_introducedY) 
summary(lmer2)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_introduced/1e+06) ~ EFN + Domatia + fixer + scale(abs_lat_native) +  
    ##     scale(total_area_native) + annual + woody + uses_num_uses +  
    ##     (1 | tribe_ncbi)
    ##    Data: legume_range_df_introducedY
    ## 
    ## REML criterion at convergence: 3717.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.3614 -0.5074  0.2251  0.6725  1.8533 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.2051   0.4529  
    ##  Residual               5.4791   2.3407  
    ## Number of obs: 814, groups:  tribe_ncbi, 39
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)              13.65018    0.40527  33.682
    ## EFN1                      0.25734    0.25026   1.028
    ## Domatia1                 -0.33832    0.97991  -0.345
    ## fixer1                   -0.51137    0.33635  -1.520
    ## scale(abs_lat_native)     0.30137    0.10549   2.857
    ## scale(total_area_native) -0.37416    0.09443  -3.962
    ## annual                    0.46486    0.30434   1.527
    ## woody                    -0.71074    0.25223  -2.818
    ## uses_num_uses             0.53144    0.04199  12.657
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   Domat1 fixer1 scl(b__) scl(t__) annual woody 
    ## EFN1        -0.049                                                     
    ## Domatia1     0.026 -0.017                                              
    ## fixer1      -0.781  0.028 -0.040                                       
    ## scl(bs_lt_)  0.015  0.064  0.061 -0.130                                
    ## scl(ttl_r_)  0.073 -0.026 -0.003 -0.109  0.060                         
    ## annual      -0.270  0.013 -0.002 -0.027 -0.017    0.183                
    ## woody       -0.552 -0.088 -0.013  0.123  0.147    0.177    0.431       
    ## uses_num_ss -0.289 -0.042 -0.019  0.089  0.035   -0.376    0.014 -0.055

``` r
Anova(lmer2, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total_area_introduced/1e+06)
    ##                              Chisq Df Pr(>Chisq)    
    ## (Intercept)              1134.4495  1  < 2.2e-16 ***
    ## EFN                         1.0574  1   0.303811    
    ## Domatia                     0.1192  1   0.729904    
    ## fixer                       2.3114  1   0.128425    
    ## scale(abs_lat_native)       8.1621  1   0.004277 ** 
    ## scale(total_area_native)   15.6999  1  7.423e-05 ***
    ## annual                      2.3330  1   0.126658    
    ## woody                       7.9404  1   0.004834 ** 
    ## uses_num_uses             160.2112  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer2)
```

![](README_files/figure-gfm/Legume%20mixed%20models-13.png)<!-- -->

``` r
#Native range size
lmer3 <- lmer(log(total_area_native/1e+6)~EFN*uses_num_uses+fixer*uses_num_uses+EFN*Domatia+fixer+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=df)
summary(lmer3)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(total_area_native/1e+06) ~ EFN * uses_num_uses + fixer *  
    ##     uses_num_uses + EFN * Domatia + fixer + annual + woody +  
    ##     scale(abs_lat_native) + uses_num_uses + (1 | tribe_ncbi)
    ##    Data: df
    ## 
    ## REML criterion at convergence: 11421.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -9.9467 -0.3532  0.1027  0.6065  2.7179 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.1127   0.3357  
    ##  Residual               1.4445   1.2019  
    ## Number of obs: 3534, groups:  tribe_ncbi, 51
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           15.04010    0.12106 124.233
    ## EFN1                   0.11520    0.10923   1.055
    ## uses_num_uses          0.18311    0.03639   5.033
    ## fixer1                -0.24345    0.10753  -2.264
    ## Domatia1               0.51836    0.30019   1.727
    ## annual                -0.07513    0.07376  -1.019
    ## woody                 -0.33316    0.05899  -5.648
    ## scale(abs_lat_native) -0.18060    0.02818  -6.410
    ## EFN1:uses_num_uses    -0.07957    0.04050  -1.965
    ## uses_num_uses:fixer1   0.12920    0.03929   3.288
    ## EFN1:Domatia1         -2.13586    0.67786  -3.151
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN1   uss_n_ fixer1 Domat1 annual woody  sc(__) EFN1:_
    ## EFN1        -0.043                                                        
    ## uses_num_ss -0.367  0.038                                                 
    ## fixer1      -0.768  0.021  0.412                                          
    ## Domatia1     0.013  0.036  0.041 -0.035                                   
    ## annual      -0.180 -0.016 -0.005 -0.032  0.011                            
    ## woody       -0.412 -0.026 -0.020  0.066 -0.009  0.392                     
    ## scl(bs_lt_)  0.035  0.020 -0.009 -0.106  0.078  0.127  0.120              
    ## EFN1:ss_nm_  0.023 -0.574 -0.097  0.012 -0.010  0.005  0.011  0.017       
    ## uss_nm_ss:1  0.352  0.005 -0.909 -0.433 -0.032  0.019 -0.012  0.015 -0.060
    ## EFN1:Domat1  0.017 -0.062 -0.017 -0.018 -0.435  0.000  0.008 -0.008 -0.052
    ##             us__:1
    ## EFN1              
    ## uses_num_ss       
    ## fixer1            
    ## Domatia1          
    ## annual            
    ## woody             
    ## scl(bs_lt_)       
    ## EFN1:ss_nm_       
    ## uss_nm_ss:1       
    ## EFN1:Domat1  0.010

``` r
Anova(lmer3, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total_area_native/1e+06)
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)           15433.8830  1  < 2.2e-16 ***
    ## EFN                       1.1124  1   0.291554    
    ## uses_num_uses            25.3275  1  4.838e-07 ***
    ## fixer                     5.1254  1   0.023578 *  
    ## Domatia                   2.9817  1   0.084210 .  
    ## annual                    1.0374  1   0.308428    
    ## woody                    31.8988  1  1.624e-08 ***
    ## scale(abs_lat_native)    41.0839  1  1.458e-10 ***
    ## EFN:uses_num_uses         3.8600  1   0.049450 *  
    ## uses_num_uses:fixer      10.8111  1   0.001009 ** 
    ## EFN:Domatia               9.9282  1   0.001628 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer3)
```

![](README_files/figure-gfm/Legume%20mixed%20models-14.png)<!-- -->

``` r
#Plot the interaction between EFN, domatia and each mutualism with human uses

#Mycorrhizae
#Successful introduction?
binomial3 <- glmer(introducedY~AM*EM+scale(total_area_native)+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(df, !is.na(myco)), family="binomial", glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(binomial3) 
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: introducedY ~ AM * EM + scale(total_area_native) + annual + woody +  
    ##     scale(abs_lat_native) + uses_num_uses + (1 | tribe_ncbi)
    ##    Data: subset(df, !is.na(myco))
    ## Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e+05))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    602.6    647.8   -291.3    582.6      671 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5318 -0.4742 -0.1826  0.4179  4.4635 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.7624   0.8732  
    ## Number of obs: 681, groups:  tribe_ncbi, 42
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.00235    0.59269  -3.378 0.000729 ***
    ## AMY                      -0.62848    0.48736  -1.290 0.197203    
    ## EMY                      -1.32063    0.76841  -1.719 0.085675 .  
    ## scale(total_area_native)  0.22286    0.14494   1.538 0.124137    
    ## annual                    0.58861    0.50982   1.155 0.248281    
    ## woody                     0.87115    0.35776   2.435 0.014892 *  
    ## scale(abs_lat_native)     0.20422    0.15065   1.356 0.175225    
    ## uses_num_uses             0.95836    0.08766  10.933  < 2e-16 ***
    ## AMY:EMY                   1.32108    0.81358   1.624 0.104423    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    scl(t__) annual woody  scl(b__) uss_n_
    ## AMY         -0.760                                                     
    ## EMY         -0.445  0.577                                              
    ## scl(ttl_r_) -0.038  0.004  0.034                                       
    ## annual      -0.270  0.089  0.041  0.162                                
    ## woody       -0.483  0.000 -0.081  0.199    0.335                       
    ## scl(bs_lt_) -0.192  0.015 -0.062  0.042    0.115  0.290                
    ## uses_num_ss -0.254 -0.061 -0.068 -0.243    0.023  0.132  0.161         
    ## AMY:EMY      0.426 -0.588 -0.876 -0.018   -0.038  0.045  0.066    0.087

``` r
Anova(binomial3, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               11.4136  1  0.0007291 ***
    ## AM                         1.6630  1  0.1972029    
    ## EM                         2.9538  1  0.0856751 .  
    ## scale(total_area_native)   2.3643  1  0.1241369    
    ## annual                     1.3330  1  0.2482811    
    ## woody                      5.9292  1  0.0148917 *  
    ## scale(abs_lat_native)      1.8377  1  0.1752252    
    ## uses_num_uses            119.5344  1  < 2.2e-16 ***
    ## AM:EM                      2.6367  1  0.1044232    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial3)
```

![](README_files/figure-gfm/Legume%20mixed%20models-15.png)<!-- -->

``` r
#No interactions with human uses

#Number of introduced ranges
lmer4 <- lmer(log(num_introduced)~EM*uses_num_uses+AM+EM+scale(total_area_native)+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(legume_range_df_introducedY, !is.na(myco)))
summary(lmer4)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(num_introduced) ~ EM * uses_num_uses + AM + EM + scale(total_area_native) +  
    ##     annual + woody + scale(abs_lat_native) + uses_num_uses +  
    ##     (1 | tribe_ncbi)
    ##    Data: subset(legume_range_df_introducedY, !is.na(myco))
    ## 
    ## REML criterion at convergence: 922.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9031 -0.7388  0.1026  0.7042  2.2842 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.08326  0.2885  
    ##  Residual               0.94283  0.9710  
    ## Number of obs: 320, groups:  tribe_ncbi, 34
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.91687    0.26669   3.438
    ## EMY                      -0.49844    0.31699  -1.572
    ## uses_num_uses             0.29792    0.02621  11.365
    ## AMY                       0.17950    0.22309   0.805
    ## scale(total_area_native) -0.17711    0.06294  -2.814
    ## annual                    0.19125    0.23369   0.818
    ## woody                    -0.35435    0.16215  -2.185
    ## scale(abs_lat_native)    -0.11853    0.06920  -1.713
    ## EMY:uses_num_uses         0.18406    0.08258   2.229
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EMY    uss_n_ AMY    scl(t__) annual woody  scl(b__)
    ## EMY         -0.194                                                     
    ## uses_num_ss -0.250  0.217                                              
    ## AMY         -0.763  0.143 -0.094                                       
    ## scl(ttl_r_) -0.038  0.049 -0.232  0.024                                
    ## annual      -0.263  0.015 -0.053  0.074  0.238                         
    ## woody       -0.443 -0.107 -0.059  0.009  0.184    0.369                
    ## scl(bs_lt_) -0.124 -0.083  0.084 -0.001  0.019    0.075  0.218         
    ## EMY:ss_nm_s  0.166 -0.733 -0.283 -0.090 -0.009   -0.011  0.013 -0.043

``` r
Anova(lmer4, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               11.8197  1  0.0005861 ***
    ## EM                         2.4725  1  0.1158531    
    ## uses_num_uses            129.1585  1  < 2.2e-16 ***
    ## AM                         0.6474  1  0.4210603    
    ## scale(total_area_native)   7.9175  1  0.0048958 ** 
    ## annual                     0.6697  1  0.4131401    
    ## woody                      4.7756  1  0.0288660 *  
    ## scale(abs_lat_native)      2.9342  1  0.0867208 .  
    ## EM:uses_num_uses           4.9671  1  0.0258337 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer4)
```

![](README_files/figure-gfm/Legume%20mixed%20models-16.png)<!-- -->

``` r
#EM uses interaction seen
#Check how many species have EM and their human uses hist

#Introduced area (if introduced at all)
lmer5 <- lmer(log(total_area_introduced/1e+6)~AM+EM+scale(total_area_native)+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(legume_range_df_introducedY, !is.na(myco)))
summary(lmer5)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_introduced/1e+06) ~ AM + EM + scale(total_area_native) +  
    ##     annual + woody + scale(abs_lat_native) + uses_num_uses +  
    ##     (1 | tribe_ncbi)
    ##    Data: subset(legume_range_df_introducedY, !is.na(myco))
    ## 
    ## REML criterion at convergence: 1467.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.0716 -0.4117  0.2491  0.6195  1.5620 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.3365   0.5801  
    ##  Residual               5.5137   2.3481  
    ## Number of obs: 320, groups:  tribe_ncbi, 34
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)              12.37901    0.62669  19.753
    ## AMY                       1.04477    0.53585   1.950
    ## EMY                       0.08356    0.51122   0.163
    ## scale(total_area_native) -0.49868    0.15079  -3.307
    ## annual                    0.99582    0.56230   1.771
    ## woody                    -0.58346    0.38550  -1.514
    ## scale(abs_lat_native)     0.43749    0.16255   2.691
    ## uses_num_uses             0.49674    0.06057   8.201
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    scl(t__) annual woody  scl(b__)
    ## AMY         -0.772                                              
    ## EMY         -0.110  0.117                                       
    ## scl(ttl_r_) -0.047  0.024  0.071                                
    ## annual      -0.274  0.074  0.015  0.233                         
    ## woody       -0.449  0.012 -0.155  0.202    0.384                
    ## scl(bs_lt_) -0.126 -0.006 -0.168  0.013    0.064  0.244         
    ## uses_num_ss -0.213 -0.125  0.022 -0.245   -0.057 -0.060  0.077

``` r
Anova(lmer5, type=3, test="F")
```

    ## Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
    ## 
    ## Response: log(total_area_introduced/1e+06)
    ##                                 F Df Df.res    Pr(>F)    
    ## (Intercept)              384.9964  1 266.73 < 2.2e-16 ***
    ## AM                         3.7822  1 305.70  0.052718 .  
    ## EM                         0.0253  1 194.01  0.873909    
    ## scale(total_area_native)  10.7231  1 300.68  0.001182 ** 
    ## annual                     3.1016  1 311.77  0.079198 .  
    ## woody                      2.2197  1 235.25  0.137600    
    ## scale(abs_lat_native)      6.8957  1 136.31  0.009628 ** 
    ## uses_num_uses             66.7292  1 309.85 7.992e-15 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer5)
```

![](README_files/figure-gfm/Legume%20mixed%20models-17.png)<!-- -->

``` r
#Native range size
lmer6 <- lmer(log(total_area_native/1e+6)~AM+EM+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(df, !is.na(myco)))
summary(lmer6)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total_area_native/1e+06) ~ AM + EM + annual + woody + scale(abs_lat_native) +  
    ##     uses_num_uses + (1 | tribe_ncbi)
    ##    Data: subset(df, !is.na(myco))
    ## 
    ## REML criterion at convergence: 2047
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.8824 -0.4290  0.1092  0.6323  2.1235 
    ## 
    ## Random effects:
    ##  Groups     Name        Variance Std.Dev.
    ##  tribe_ncbi (Intercept) 0.09374  0.3062  
    ##  Residual               1.10887  1.0530  
    ## Number of obs: 681, groups:  tribe_ncbi, 42
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           15.59693    0.19010  82.046
    ## AMY                    0.06331    0.15018   0.422
    ## EMY                    0.14946    0.14122   1.058
    ## annual                -0.67965    0.18476  -3.679
    ## woody                 -0.61038    0.13021  -4.688
    ## scale(abs_lat_native) -0.02858    0.05494  -0.520
    ## uses_num_uses          0.16579    0.01909   8.685
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) AMY    EMY    annual woody  sc(__)
    ## AMY         -0.728                                   
    ## EMY         -0.219  0.234                            
    ## annual      -0.269  0.037  0.012                     
    ## woody       -0.519  0.005 -0.123  0.344              
    ## scl(bs_lt_) -0.167  0.010 -0.038  0.097  0.260       
    ## uses_num_ss -0.166 -0.073  0.038  0.020  0.006  0.096

``` r
Anova(lmer6, type=3, test="F")
```

    ## Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)
    ## 
    ## Response: log(total_area_native/1e+06)
    ##                               F Df Df.res    Pr(>F)    
    ## (Intercept)           6635.7694  1 324.27 < 2.2e-16 ***
    ## AM                       0.1763  1 671.16 0.6747332    
    ## EM                       1.0931  1 472.28 0.2963181    
    ## annual                  13.4563  1 673.74 0.0002635 ***
    ## woody                   21.4553  1 403.51 4.893e-06 ***
    ## scale(abs_lat_native)    0.2630  1 272.10 0.6084558    
    ## uses_num_uses           75.0464  1 673.62 < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer6)
```

![](README_files/figure-gfm/Legume%20mixed%20models-18.png)<!-- -->

#### Multiple mutualisms

We can also look at the effects of multiple mutualisms.

``` r
df$AM.num <- ifelse(df$AM == "Y", 1, ifelse(df$AM == "N", 0, NA))
df$EM.num <- ifelse(df$EM == "Y", 1, ifelse(df$EM == "N", 0, NA))

df$num_mutualisms <- as.numeric(as.character(df$fixer))+as.numeric(as.character(df$EFN))+as.numeric(as.character(df$AM.num))+as.numeric(as.character(df$Domatia))+as.numeric(as.character(df$EM.num))

summary.mnum2 <- ungroup(subset(df, !is.na(num_mutualisms)) %>% group_by(num_mutualisms, introducedY) %>% dplyr::summarize(n=n()))
summary.mnum2.wide <- spread(summary.mnum2, key = introducedY, value=n)
colnames(summary.mnum2.wide) <- c("num_mutualisms","Not_introduced",  "Introduced")
summary.mnum2.wide$total <- summary.mnum2.wide$Not_introduced+summary.mnum2.wide$Introduced
summary.mnum2.wide$prop.introduced <- summary.mnum2.wide$Introduced/(summary.mnum2.wide$total)
prop.mnum <- paste0(summary.mnum2.wide$Introduced, "/", summary.mnum2.wide$total)

inset_p_mnum <- ggplot(data=subset(summary.mnum2.wide, !is.na(num_mutualisms)), aes(x=num_mutualisms, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+xlab("Mutualisms (no.)")+ylab("Introduced (prop.)")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 2, y = 0.55, label = "**")+geom_text(aes(x=num_mutualisms, y=0.05, label=prop.mnum), color="white")

summary.mnum <- subset(df, !is.na(num_mutualisms) & introducedY == 1) %>% group_by(num_mutualisms) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n))
summary.mnum$num_mutualisms <- as.factor(summary.mnum$num_mutualisms)

p_num_mutualisms <- ggplot(data=summary.mnum, aes(x=num_mutualisms, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num_mutualisms, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+theme_cowplot()+ylab("Introduced ranges (no.)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num_mutualisms, y= 3, label=n))+annotate("text", x=3, y=15.5, label="p = 0.06")
p_num_mutualisms
```

![](README_files/figure-gfm/Number%20of%20mutualisms-1.png)<!-- -->

``` r
lmer7 <-  lmer(log(num_introduced)~num_mutualisms+annual+scale(abs_lat_native)+woody+uses_num_uses+scale(total_area_native)+(1|genus), data=subset(df, !is.na(num_mutualisms) & introducedY == 1))
summary(lmer7)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(num_introduced) ~ num_mutualisms + annual + scale(abs_lat_native) +  
    ##     woody + uses_num_uses + scale(total_area_native) + (1 | genus)
    ##    Data: subset(df, !is.na(num_mutualisms) & introducedY == 1)
    ## 
    ## REML criterion at convergence: 928.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7814 -0.7254  0.1240  0.6769  2.2338 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  genus    (Intercept) 0.08196  0.2863  
    ##  Residual             0.94288  0.9710  
    ## Number of obs: 320, groups:  genus, 135
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.70320    0.24419   2.880
    ## num_mutualisms            0.15926    0.09372   1.699
    ## annual                    0.21943    0.23495   0.934
    ## scale(abs_lat_native)    -0.08893    0.06355  -1.399
    ## woody                    -0.37071    0.15853  -2.338
    ## uses_num_uses             0.32070    0.02532  12.664
    ## scale(total_area_native) -0.17052    0.06279  -2.716
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_mtl annual scl(b__) woody  uss_n_
    ## num_mutlsms -0.765                                     
    ## annual      -0.262  0.026                              
    ## scl(bs_lt_) -0.115 -0.045  0.051                       
    ## woody       -0.450 -0.015  0.405  0.271                
    ## uses_num_ss -0.304 -0.035 -0.043  0.090   -0.059       
    ## scl(ttl_r_) -0.028 -0.026  0.219  0.015    0.247 -0.243

``` r
Anova(lmer7, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(num_introduced)
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)                8.2926  1   0.003981 ** 
    ## num_mutualisms             2.8880  1   0.089241 .  
    ## annual                     0.8722  1   0.350344    
    ## scale(abs_lat_native)      1.9581  1   0.161722    
    ## woody                      5.4682  1   0.019366 *  
    ## uses_num_uses            160.3828  1  < 2.2e-16 ***
    ## scale(total_area_native)   7.3743  1   0.006616 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer7)
```

![](README_files/figure-gfm/Number%20of%20mutualisms-2.png)<!-- -->

``` r
binomial4 <-  glmer(introducedY~num_mutualisms+annual+woody+scale(abs_lat_native)+uses_num_uses+scale(total_area_native)+(1|genus), data=subset(df, !is.na(num_mutualisms)), family="binomial")
summary(binomial4)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ num_mutualisms + annual + woody + scale(abs_lat_native) +  
    ##     uses_num_uses + scale(total_area_native) + (1 | genus)
    ##    Data: subset(df, !is.na(num_mutualisms))
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    616.8    653.0   -300.4    600.8      673 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6612 -0.4410 -0.2309  0.4093  4.1742 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  genus  (Intercept) 0.7663   0.8754  
    ## Number of obs: 681, groups:  genus, 235
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -3.40743    0.53923  -6.319 2.63e-10 ***
    ## num_mutualisms            0.50841    0.19826   2.564   0.0103 *  
    ## annual                    0.92728    0.50582   1.833   0.0668 .  
    ## woody                     0.62268    0.35660   1.746   0.0808 .  
    ## scale(abs_lat_native)     0.13743    0.13893   0.989   0.3226    
    ## uses_num_uses             0.99133    0.09623  10.302  < 2e-16 ***
    ## scale(total_area_native)  0.28482    0.14761   1.930   0.0537 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) nm_mtl annual woody  scl(b__) uss_n_
    ## num_mutlsms -0.755                                     
    ## annual      -0.278  0.024                              
    ## woody       -0.523 -0.002  0.385                       
    ## scl(bs_lt_) -0.187 -0.056  0.073  0.374                
    ## uses_num_ss -0.432  0.110  0.071  0.125  0.194         
    ## scl(ttl_r_) -0.027 -0.057  0.160  0.261  0.025   -0.220
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## Model failed to converge with max|grad| = 0.0066332 (tol = 0.002, component 1)

``` r
Anova(binomial4, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               39.9308  1  2.631e-10 ***
    ## num_mutualisms             6.5759  1    0.01034 *  
    ## annual                     3.3606  1    0.06677 .  
    ## woody                      3.0490  1    0.08079 .  
    ## scale(abs_lat_native)      0.9784  1    0.32259    
    ## uses_num_uses            106.1225  1  < 2.2e-16 ***
    ## scale(total_area_native)   3.7233  1    0.05366 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial4)
```

![](README_files/figure-gfm/Number%20of%20mutualisms-3.png)<!-- -->

``` r
#Rerun with optimization

fig3 <- plot_grid(inset_p_mnum, p_num_mutualisms, ncol=1, labels="AUTO")
fig3
```

![](README_files/figure-gfm/Number%20of%20mutualisms-4.png)<!-- -->

``` r
save_plot("Figure3.pdf", fig3, base_height = 8, base_width=5)

sum.nmum <- subset(df, !is.na(num_mutualisms)) %>% group_by(AM, EM, Domatia, EFN, fixer, num_mutualisms) %>% dplyr::summarize(n=n())
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
ncontig <- read.csv("Ant_noncontig.csv") #number of non-contiguous introduced ranges for INTRODUCED ants only
abslat <- read.csv("absolute_native_lat_ants7Feb.csv") #Absolute midpoint latitude of native range

#merging invaded and native area datasets
nat_inv_area <- merge(invarea, antarea, by='Phy', all.y=TRUE) 
nat_inv_area <- merge(nat_inv_area, ncontig, by='Phy', all.x=TRUE)
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
colnames(area)[[15]] <- "ExoticY"
area <- merge(area, indoor[, c("Phy", "Exotic_status")], by="Phy", all=TRUE)
colnames(area)[[16]] <- "IndoorY"
area <- merge(area, intercepted[, c("Phy", "Exotic_status")], by="Phy", all=TRUE)
colnames(area)[[17]] <- "InterceptedY"

#Introduced? 
area$introducedY <- ifelse(area$total.area.introduced > 0, 1, 0)
area$IndoorY <- ifelse(!is.na(area$IndoorY) & area$IndoorY == "Indoor Introduced" & area$introducedY == 1, 1, ifelse(is.na(area$IndoorY) & area$introducedY == 1, 0, NA))
area$ExoticY <- ifelse(!is.na(area$ExoticY) & area$ExoticY == "Exotic" & area$introducedY == 1, 1, ifelse(is.na(area$ExoticY) & area$introducedY == 1, 0, NA))
area$InterceptedY <- ifelse(!is.na(area$InterceptedY) & area$InterceptedY == "Intercepted" & area$introducedY == 1, 1, ifelse(is.na(area$InterceptedY) & area$introducedY == 1, 0, NA))
```

We will use the ant subfamily as a random effect in models, and so need
to download it for each ant genus in the dataset.

``` r
#Not run
#Code included for reproducibility
ant <- data.frame(genus = unique(word(area$Phy, 1, 1, sep="_"))) #creating the vector of unique genera
ant$subfamily <- NA #create an empty vector ant subfamily 
#Loop through all ant genera
for (i in 1:length(ant$genus)){
  ant[i, "subfamily"] <- tax_name(sci= genus[i], get = 'subfamily', db = 'itis')$subfamily
}
write.csv(ant, "subfamily.csv") 

#Merge ant subfamily with ant range dataset
area$genus <- word(area$Phy, 1, 1, sep="_")
area <- merge(area, ant, by="genus", all.x=TRUE)
write.csv(area, "ant_areas.csv")
```

### Summarize ant dataset

``` r
area <- read.csv("ant_areas.csv")
area$subfamily <- as.factor(area$subfamily)

sum(complete.cases(subset(area, !is.na(introducedY))[, c("EFN", "Domatia", "Seed_Dispersal")])) #Total taxa with at least some trait data 
```

    ## [1] 3023

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(EFN) %>% dplyr::summarize(n=n()))
```

| EFN |     n |
|----:|------:|
|   0 |  2910 |
|   1 |   115 |
|  NA | 12271 |

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(Domatia) %>% dplyr::summarize(n=n()))
```

| Domatia |     n |
|--------:|------:|
|       0 |  2967 |
|       1 |    58 |
|      NA | 12271 |

``` r
kable(subset(area, !is.na(introducedY)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n()))
```

| Seed_Dispersal |     n |
|---------------:|------:|
|              0 |  2758 |
|              1 |   297 |
|             NA | 12241 |

``` r
pt_size <- 3
y_limits <- c(-100000, 1.5e+7)
er_width <- 0.1
y_text <- -100000

summary.ant.efn <- subset(area, introducedY == 1 & !is.na(EFN)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

summary.ant.efn1 <- subset(area, introducedY == 1 & !is.na(EFN)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antEFN <- ggplot(data=summary.ant.efn, aes(x=EFN, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced range (sq. km)")+xlab("Visits EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "*")

p_antEFN1 <- ggplot(data=summary.ant.efn1, aes(x=EFN, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Number of non-contiguous introduced ranges")+xlab("Visits EFNs")#+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 5, label = "*")
#Adjust labels

summary.ant.efn2 <- ungroup(subset(area, !is.na(EFN)) %>% group_by(EFN, introducedY) %>% dplyr::summarize(n=n()))
```

    ## `summarise()` has grouped output by 'EFN'. You can override using the `.groups`
    ## argument.

``` r
summary.ant.efn2.wide <- spread(summary.ant.efn2, key = introducedY, value=n)
colnames(summary.ant.efn2.wide) <- c("EFN","Not_introduced",  "Introduced")
summary.ant.efn2.wide$total <- summary.ant.efn2.wide$Not_introduced+summary.ant.efn2.wide$Introduced
summary.ant.efn2.wide$prop.introduced <- summary.ant.efn2.wide$Introduced/(summary.ant.efn2.wide$total)
prop.ant.efn <- paste0(summary.ant.efn2.wide$Introduced, "/", summary.ant.efn2.wide$total)

inset_p_antEFN <- ggplot(data=summary.ant.efn2.wide, aes(x=EFN, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Visits EFNs")#+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=EFN, y=0.05, label=prop.ant.efn), color="white")

summary.ant.dom <- subset(area, introducedY == 1 & !is.na(Domatia)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

summary.ant.dom1 <- subset(area, introducedY == 1 & !is.na(Domatia)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antdom <- ggplot(data=summary.ant.dom1, aes(x=Domatia, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("No. of non-contiguous ranges")+xlab("Nests in domatia")#+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "ns")

summary.ant.dom2 <- ungroup(subset(area, !is.na(Domatia)) %>% group_by(Domatia, introducedY) %>% dplyr::summarize(n=n()))
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

summary.ant.elaiosome <- subset(area, introducedY == 1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n(), mean_introduced = mean(total.area.introduced, na.rm=TRUE), sd_introduced = sd(total.area.introduced, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

summary.ant.elaiosome1 <- subset(area, introducedY == 1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_elaiosome <- ggplot(data=summary.ant.elaiosome1, aes(x=Seed_Dispersal, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Seed_Dispersal, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("No. of non-contiguous ranges")+xlab("Disperses seeds")+geom_text(aes(x=Seed_Dispersal, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))#+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7e+6, label = "ns")

summary.ant.seed2 <- ungroup(subset(area, !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, introducedY) %>% dplyr::summarize(n=n()))
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
fig4bottom <- plot_grid(p_antEFN1, p_elaiosome, p_antdom, nrow=1, labels=c("D", "E", "F"))
fig4 <- plot_grid(fig4top, fig4bottom, nrow=2)
fig4
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
save_plot("Figure4i.pdf", fig4, base_width=8, base_height =8)
```

#### Native range size

``` r
y_limits = c(3200000, 1.25e+07)
y_text = 3300000
summary.ant.efn.native <- subset(area, !is.na(EFN) & !is.na(total.area.native)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antEFN_native <- ggplot(data=summary.ant.efn.native, aes(x=EFN, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Visits EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 8100000, label = "***")

summary.ant.dom.native <- subset(area, !is.na(Domatia)& !is.na(total.area.native)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antdom_native <- ggplot(data=summary.ant.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nests in domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 8100000, label = "*")

summary.ant.seed.native <- subset(area, !is.na(Seed_Dispersal)& !is.na(total.area.native)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

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
area$num.mm <- as.numeric(as.character(area$Seed_Dispersal))+as.numeric(as.character(area$Domatia))+as.numeric(as.character(area$EFN))

binomial10 <- glmer(introducedY ~ num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=area, family ="binomial")
summary(binomial10)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ num.mm + scale(abs_lat_native) + scale(total.area.native) +  
    ##     (1 | subfamily)
    ##    Data: area
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1762.9   1792.7   -876.4   1752.9     2841 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.1257 -0.3286 -0.2687 -0.2193  5.4137 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.4536   0.6735  
    ## Number of obs: 2846, groups:  subfamily, 16
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -3.04409    0.33296  -9.142  < 2e-16 ***
    ## num.mm                    0.81645    0.12290   6.643 3.08e-11 ***
    ## scale(abs_lat_native)    -0.25739    0.05660  -4.547 5.43e-06 ***
    ## scale(total.area.native)  0.36587    0.03318  11.026  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) num.mm sc(__)
    ## num.mm      -0.027              
    ## scl(bs_lt_)  0.113  0.035       
    ## scl(ttl.r.) -0.140 -0.260 -0.267

``` r
Anova(binomial10, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)               83.585  1  < 2.2e-16 ***
    ## num.mm                    44.128  1  3.075e-11 ***
    ## scale(abs_lat_native)     20.678  1  5.433e-06 ***
    ## scale(total.area.native) 121.562  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial10)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
lmer11 <- lmer(log(n_introduced_ranges)~num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1))
summary(lmer11)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(n_introduced_ranges) ~ num.mm + scale(abs_lat_native) + scale(total.area.native) +  
    ##     (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ## REML criterion at convergence: 1055.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5344 -0.7106 -0.2446  0.4045  2.7743 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.03355  0.1832  
    ##  Residual              1.47881  1.2161  
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.87702    0.12209   7.183
    ## num.mm                    0.21928    0.11158   1.965
    ## scale(abs_lat_native)    -0.35900    0.06682  -5.373
    ## scale(total.area.native)  0.12184    0.06751   1.805
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) num.mm sc(__)
    ## num.mm      -0.330              
    ## scl(bs_lt_)  0.030 -0.022       
    ## scl(ttl.r.)  0.029 -0.375 -0.132

``` r
Anova(lmer11, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(n_introduced_ranges)
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)              51.5973  1  6.813e-13 ***
    ## num.mm                    3.8624  1    0.04938 *  
    ## scale(abs_lat_native)    28.8646  1  7.762e-08 ***
    ## scale(total.area.native)  3.2569  1    0.07113 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary.ant.mm <- ungroup(subset(area, !is.na(introducedY)) %>% group_by(num.mm, introducedY) %>% dplyr::summarize(n=n()))
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

summary.ant.mnum <- subset(area, !is.na(num.mm) & introducedY == 1) %>% group_by(num.mm) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(total.area.introduced, na.rm=TRUE), sd_area_introduced = sd(total.area.introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n))
summary.ant.mnum$num.mm <- as.factor(summary.ant.mnum$num.mm)

summary.ant.mnum1 <- subset(area, !is.na(num.mm) & introducedY == 1) %>% group_by(num.mm) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_area_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n))
summary.ant.mnum$num.mm <- as.factor(summary.ant.mnum$num.mm)


p_num_mm <- ggplot(data=summary.ant.mnum1, aes(x=num.mm, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num.mm, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+theme_cowplot()+ylab("Introduced range (sq. km)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num.mm, y= -1, label=n))+annotate("text", x=2.5, y=10, label="p = 0.06")

fig6 <- plot_grid(inset_p_antmm, p_num_mm, ncol=1, labels="AUTO")
fig6
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
save_plot("Figure6.pdf", fig6, base_height = 8, base_width=5)
```

### Mixed models

``` r
#Correlations between variables
df_ant_cor <- data.frame(area$abs_lat_native, as.numeric(area$EFN), as.numeric(area$Domatia), as.numeric(area$Seed_Dispersal))
  #covariates
# <- cor(df_cor)
colnames(df_ant_cor) <- c("Abs_lat", "EFN", "Domatia", "Seed_dispersal")
corr_ant <- rcorr(as.matrix(df_ant_cor))

plotr <- corrplot(corr_ant$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) #correlation coeffs
```

![](README_files/figure-gfm/ant%20glmms-1.png)<!-- -->

``` r
plotp <- corrplot(corr_ant$P, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45) #p-values
```

![](README_files/figure-gfm/ant%20glmms-2.png)<!-- -->

``` r
#Introduction success
binomial6 <- lme4::glmer(introducedY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=area, family="binomial")
summary(binomial6)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## introducedY ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | subfamily)
    ##    Data: area
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   1751.9   1793.6   -868.9   1737.9     2839 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.2153 -0.3266 -0.2677 -0.2157  5.6050 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.406    0.6372  
    ## Number of obs: 2846, groups:  subfamily, 16
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.97877    0.33203  -8.972  < 2e-16 ***
    ## EFN                       0.95525    0.24452   3.907 9.36e-05 ***
    ## Domatia                  -0.65227    0.46376  -1.406     0.16    
    ## Seed_Dispersal            1.07159    0.17973   5.962 2.49e-09 ***
    ## scale(abs_lat_native)    -0.29069    0.05826  -4.990 6.04e-07 ***
    ## scale(total.area.native)  0.35759    0.03350  10.674  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Domati Sd_Dsp sc(__)
    ## EFN          0.011                            
    ## Domatia     -0.020 -0.134                     
    ## Seed_Dsprsl -0.001 -0.142 -0.001              
    ## scl(bs_lt_)  0.109  0.106  0.098 -0.144       
    ## scl(ttl.r.) -0.144 -0.102 -0.046 -0.217 -0.236

``` r
car::Anova(binomial6, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: introducedY
    ##                             Chisq Df Pr(>Chisq)    
    ## (Intercept)               80.4880  1  < 2.2e-16 ***
    ## EFN                       15.2622  1  9.357e-05 ***
    ## Domatia                    1.9782  1     0.1596    
    ## Seed_Dispersal            35.5466  1  2.490e-09 ***
    ## scale(abs_lat_native)     24.8995  1  6.040e-07 ***
    ## scale(total.area.native) 113.9236  1  < 2.2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial6)
```

![](README_files/figure-gfm/ant%20glmms-3.png)<!-- -->

``` r
#Update optimizer

#Total introduced area
lmer9 <- lmer(log(total.area.introduced)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1))
summary(lmer9)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: log(total.area.introduced) ~ EFN + Domatia + Seed_Dispersal +  
    ##     scale(abs_lat_native) + scale(total.area.native) + (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ## REML criterion at convergence: 1459.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -4.2544 -0.6217  0.0733  0.6040  2.2167 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.1413   0.376   
    ##  Residual              5.3134   2.305   
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               12.5140     0.2415  51.823
    ## EFN                        0.9154     0.4132   2.216
    ## Domatia                   -0.1931     0.9200  -0.210
    ## Seed_Dispersal            -0.1241     0.3316  -0.374
    ## scale(abs_lat_native)     -0.3618     0.1323  -2.736
    ## scale(total.area.native)   0.2280     0.1291   1.766
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Domati Sd_Dsp sc(__)
    ## EFN         -0.126                            
    ## Domatia     -0.049 -0.232                     
    ## Seed_Dsprsl -0.246 -0.228  0.038              
    ## scl(bs_lt_)  0.047  0.146  0.115 -0.223       
    ## scl(ttl.r.)  0.035 -0.120 -0.036 -0.330 -0.089

``` r
Anova(lmer9, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total.area.introduced)
    ##                              Chisq Df Pr(>Chisq)    
    ## (Intercept)              2685.6120  1    < 2e-16 ***
    ## EFN                         4.9087  1    0.02672 *  
    ## Domatia                     0.0441  1    0.83375    
    ## Seed_Dispersal              0.1400  1    0.70832    
    ## scale(abs_lat_native)       7.4854  1    0.00622 ** 
    ## scale(total.area.native)    3.1194  1    0.07737 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer9)
```

![](README_files/figure-gfm/ant%20glmms-4.png)<!-- -->

``` r
#Total non-contiguous ranges
lmer9i <- lmer(log(n_introduced_ranges)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1))
summary(lmer9i)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(n_introduced_ranges) ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ## REML criterion at convergence: 1048.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8157 -0.6595 -0.2706  0.3880  2.6915 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.06417  0.2533  
    ##  Residual              1.44778  1.2032  
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error t value
    ## (Intercept)               0.86834    0.14557   5.965
    ## EFN                       0.70595    0.21620   3.265
    ## Domatia                   0.03495    0.48059   0.073
    ## Seed_Dispersal           -0.13332    0.17344  -0.769
    ## scale(abs_lat_native)    -0.31178    0.06930  -4.499
    ## scale(total.area.native)  0.14733    0.06758   2.180
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Domati Sd_Dsp sc(__)
    ## EFN         -0.106                            
    ## Domatia     -0.047 -0.231                     
    ## Seed_Dsprsl -0.206 -0.229  0.038              
    ## scl(bs_lt_)  0.048  0.150  0.114 -0.223       
    ## scl(ttl.r.)  0.027 -0.116 -0.035 -0.332 -0.086

``` r
Anova(lmer9i, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(n_introduced_ranges)
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)              35.5841  1  2.443e-09 ***
    ## EFN                      10.6618  1   0.001094 ** 
    ## Domatia                   0.0053  1   0.942031    
    ## Seed_Dispersal            0.5909  1   0.442080    
    ## scale(abs_lat_native)    20.2410  1  6.827e-06 ***
    ## scale(total.area.native)  4.7523  1   0.029259 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer9i)
```

![](README_files/figure-gfm/ant%20glmms-5.png)<!-- -->

``` r
#Native range area
lmer10 <- lmer(log(total.area.native)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+(1|subfamily), data=area)
summary(lmer10)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: 
    ## log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     (1 | subfamily)
    ##    Data: area
    ## 
    ## REML criterion at convergence: 10911.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -6.5375 -0.4862  0.1723  0.6810  2.0523 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.8845   0.9405  
    ##  Residual              2.6646   1.6323  
    ## Number of obs: 2846, groups:  subfamily, 16
    ## 
    ## Fixed effects:
    ##                       Estimate Std. Error t value
    ## (Intercept)           14.23236    0.27307  52.119
    ## EFN                    0.67473    0.16313   4.136
    ## Domatia                0.45725    0.23003   1.988
    ## Seed_Dispersal         1.34914    0.11386  11.849
    ## scale(abs_lat_native)  0.13321    0.02811   4.739
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Domati Sd_Dsp
    ## EFN         -0.008                     
    ## Domatia     -0.024 -0.119              
    ## Seed_Dsprsl -0.014 -0.202  0.003       
    ## scl(bs_lt_)  0.005  0.068  0.100 -0.142

``` r
Anova(lmer10, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: log(total.area.native)
    ##                           Chisq Df Pr(>Chisq)    
    ## (Intercept)           2716.4068  1  < 2.2e-16 ***
    ## EFN                     17.1083  1  3.531e-05 ***
    ## Domatia                  3.9513  1    0.04684 *  
    ## Seed_Dispersal         140.3928  1  < 2.2e-16 ***
    ## scale(abs_lat_native)   22.4626  1  2.143e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(lmer10)
```

![](README_files/figure-gfm/ant%20glmms-6.png)<!-- -->

### By introduction mode

``` r
# Exotic
area$ExoticY <- as.factor(area$ExoticY)
exotic.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, ExoticY) %>% dplyr::summarize(n=n()))
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

exotic.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, ExoticY) %>% dplyr::summarize(n=n()))
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

exotic.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, ExoticY) %>% dplyr::summarize(n=n()))
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

binomial7 <- glmer(ExoticY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1), family="binomial")
summary(binomial7)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: ExoticY ~ EFN + Domatia + Seed_Dispersal + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    426.0    452.5   -206.0    412.0      316 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5785 -1.0221  0.5983  0.7440  1.5618 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.3139   0.5603  
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)               0.51199    0.30234   1.693   0.0904 .
    ## EFN                       0.66664    0.39923   1.670   0.0950 .
    ## Domatia                  -0.70096    0.83700  -0.837   0.4023  
    ## Seed_Dispersal           -0.29867    0.30682  -0.973   0.3303  
    ## scale(abs_lat_native)    -0.18516    0.12365  -1.497   0.1343  
    ## scale(total.area.native)  0.08487    0.11997   0.707   0.4793  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Domati Sd_Dsp sc(__)
    ## EFN         -0.063                            
    ## Domatia     -0.047 -0.242                     
    ## Seed_Dsprsl -0.188 -0.220  0.039              
    ## scl(bs_lt_)  0.056  0.155  0.119 -0.224       
    ## scl(ttl.r.)  0.044 -0.108 -0.035 -0.345 -0.082

``` r
Anova(binomial7, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: ExoticY
    ##                           Chisq Df Pr(>Chisq)  
    ## (Intercept)              2.8676  1    0.09038 .
    ## EFN                      2.7883  1    0.09495 .
    ## Domatia                  0.7014  1    0.40233  
    ## Seed_Dispersal           0.9476  1    0.33033  
    ## scale(abs_lat_native)    2.2424  1    0.13427  
    ## scale(total.area.native) 0.5005  1    0.47929  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial7)
```

![](README_files/figure-gfm/Introduction%20mode-1.png)<!-- -->

``` r
#Indoors?
area$IndoorY <- as.factor(area$IndoorY)

indoor.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, IndoorY) %>% dplyr::summarize(n=n()))
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

indoor.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, IndoorY) %>% dplyr::summarize(n=n()))
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

indoor.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, IndoorY) %>% dplyr::summarize(n=n()))
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

binomial8 <- glmer(IndoorY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1), family="binomial")
```

    ## boundary (singular) fit: see help('isSingular')

``` r
summary(binomial8)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: IndoorY ~ EFN + Seed_Dispersal + Domatia + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    373.9    400.3   -179.9    359.9      316 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9927 -0.6204 -0.5154  0.9101  2.9157 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0        0       
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               -0.9708     0.1516  -6.405  1.5e-10 ***
    ## EFN                        1.1327     0.3810   2.973  0.00295 ** 
    ## Seed_Dispersal            -0.9597     0.3704  -2.591  0.00958 ** 
    ## Domatia                   -0.2207     0.8283  -0.266  0.78994    
    ## scale(abs_lat_native)     -0.1760     0.1324  -1.329  0.18383    
    ## scale(total.area.native)   0.3188     0.1244   2.562  0.01041 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Sd_Dsp Domati sc(__)
    ## EFN         -0.206                            
    ## Seed_Dsprsl -0.334 -0.308                     
    ## Domatia     -0.052 -0.241  0.043              
    ## scl(bs_lt_)  0.119  0.106 -0.175  0.125       
    ## scl(ttl.r.)  0.005 -0.076 -0.369 -0.038 -0.144
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
Anova(binomial8, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: IndoorY
    ##                            Chisq Df Pr(>Chisq)    
    ## (Intercept)              41.0230  1  1.504e-10 ***
    ## EFN                       8.8408  1   0.002946 ** 
    ## Seed_Dispersal            6.7115  1   0.009579 ** 
    ## Domatia                   0.0710  1   0.789943    
    ## scale(abs_lat_native)     1.7664  1   0.183834    
    ## scale(total.area.native)  6.5629  1   0.010412 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(binomial8)
```

![](README_files/figure-gfm/Introduction%20mode-2.png)<!-- -->

``` r
#Intercepted
area$InterceptedY <- as.factor(area$InterceptedY)

intercepted.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, InterceptedY) %>% dplyr::summarize(n=n()))
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

intercepted.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, InterceptedY) %>% dplyr::summarize(n=n()))
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

intercepted.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, InterceptedY) %>% dplyr::summarize(n=n()))
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

binomial9 <- glmer(InterceptedY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|subfamily), data=subset(area, introducedY == 1), family="binomial")
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model is nearly unidentifiable: large eigenvalue ratio
    ##  - Rescale variables?

``` r
summary(binomial9)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: 
    ## InterceptedY ~ EFN + Seed_Dispersal + Domatia + scale(abs_lat_native) +  
    ##     scale(total.area.native) + (1 | subfamily)
    ##    Data: subset(area, introducedY == 1)
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    410.3    436.7   -198.1    396.3      316 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9205 -0.9088  0.4567  0.8354  2.2000 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  subfamily (Intercept) 0.3681   0.6067  
    ## Number of obs: 323, groups:  subfamily, 8
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)               -0.0409     0.3145  -0.130  0.89654   
    ## EFN                        0.4911     0.4507   1.090  0.27584   
    ## Seed_Dispersal             0.9842     0.3336   2.951  0.00317 **
    ## Domatia                   16.3172   512.0006   0.032  0.97458   
    ## scale(abs_lat_native)     -0.4317     0.1337  -3.230  0.00124 **
    ## scale(total.area.native)   0.1212     0.1328   0.913  0.36112   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) EFN    Sd_Dsp Domati sc(__)
    ## EFN         -0.067                            
    ## Seed_Dsprsl -0.164 -0.209                     
    ## Domatia      0.000  0.001 -0.001              
    ## scl(bs_lt_)  0.065  0.153 -0.305  0.001       
    ## scl(ttl.r.)  0.025 -0.072 -0.303  0.000 -0.053
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## Model is nearly unidentifiable: large eigenvalue ratio
    ##  - Rescale variables?

``` r
Anova(binomial9, type=3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: InterceptedY
    ##                            Chisq Df Pr(>Chisq)   
    ## (Intercept)               0.0169  1   0.896541   
    ## EFN                       1.1875  1   0.275836   
    ## Seed_Dispersal            8.7055  1   0.003173 **
    ## Domatia                   0.0010  1   0.974576   
    ## scale(abs_lat_native)    10.4298  1   0.001240 **
    ## scale(total.area.native)  0.8340  1   0.361119   
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

area$n_introduced_ranges <- ifelse(area$total.area.introduced == 0, 0, area$n_introduced_ranges)

area_phy <- area[area$Phy %in% phy_int, c("Phy", "abs_lat_native", "total.area.introduced", "total.area.native", "n_introduced_ranges", "EFN", "Domatia", "Seed_Dispersal")] #dataset for pgls
area_phy <- area_phy[complete.cases(area_phy), ]
#write.csv(efn_area_phy, "Ant_list_EFN.csv")
#efn_area_phy <- efn_area_phy[match(pruned_ant_tree$tip.label, efn_area_phy$Phy),]

#PGLS model
#a1 <- gls(log(total.area.native) ~ EFN + Domatia + Seed_Dispersal + EFN*Domatia + Domatia*Seed_Dispersal + Seed_Dispersal*EFN + abs_lat_native, 
                #correlation = corPagel(1, phy = pruned_ant_tree, form = ~ Phy), 
                #method = "ML", data = area_phy) 
#summary(a1) 

#removing interactions
a11 <- gls(((total.area.native)^(1/4)) ~ EFN + Domatia + Seed_Dispersal + abs_lat_native, 
                correlation = corPagel(1, phy = pruned_ant_tree, form = ~ Phy), 
                method = "ML", data = area_phy) 
summary(a11)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: ((total.area.native)^(1/4)) ~ EFN + Domatia + Seed_Dispersal +      abs_lat_native 
    ##   Data: area_phy 
    ##        AIC      BIC    logLik
    ##   6064.151 6096.666 -3025.075
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.4245393 
    ## 
    ## Coefficients:
    ##                   Value Std.Error   t-value p-value
    ## (Intercept)    40.11204  3.403817 11.784429  0.0000
    ## EFN             7.14600  1.667334  4.285886  0.0000
    ## Domatia        -1.24239  2.452935 -0.506492  0.6127
    ## Seed_Dispersal 10.35957  1.252365  8.272005  0.0000
    ## abs_lat_native  0.03896  0.040162  0.969951  0.3324
    ## 
    ##  Correlation: 
    ##                (Intr) EFN    Domati Sd_Dsp
    ## EFN            -0.030                     
    ## Domatia        -0.036 -0.142              
    ## Seed_Dispersal -0.018 -0.192  0.034       
    ## abs_lat_native -0.189  0.032  0.078 -0.036
    ## 
    ## Standardized residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.4956861 -0.4287396  0.2710189  0.8615139  2.8378220 
    ## 
    ## Residual standard error: 14.2025 
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
    ## (Intercept)        1 158.43851  <.0001
    ## EFN                1  34.70975  <.0001
    ## Domatia            1   0.79250  0.3736
    ## Seed_Dispersal     1  69.08979  <.0001
    ## abs_lat_native     1   0.94080  0.3324

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
    ## EFN                1.9445287 0.6424026  3.026963  0.0026
    ## Domatia           -0.6669233 0.9330685 -0.714764  0.4750
    ## Seed_Dispersal     2.4630599 0.4944578  4.981335  0.0000
    ## total.area.native  0.0000002 0.0000000  7.260721  0.0000
    ## abs_lat_native    -0.0693057 0.0155102 -4.468397  0.0000
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN    Domati Sd_Dsp ttl.r.
    ## EFN               -0.013                            
    ## Domatia           -0.032 -0.141                     
    ## Seed_Dispersal     0.005 -0.121  0.025              
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

``` r
#Non-contiguous range
a3 <- gls(log(n_introduced_ranges+1) ~  EFN + Domatia + Seed_Dispersal +  total.area.native + abs_lat_native,  correlation = corPagel(1, phy = pruned_ant_tree, form = ~Phy), method = "ML", data = area_phy)  
summary(a3)
```

    ## Generalized least squares fit by maximum likelihood
    ##   Model: log(n_introduced_ranges + 1) ~ EFN + Domatia + Seed_Dispersal +      total.area.native + abs_lat_native 
    ##   Data: area_phy 
    ##        AIC      BIC    logLik
    ##   1806.334 1843.494 -895.1668
    ## 
    ## Correlation Structure: corPagel
    ##  Formula: ~Phy 
    ##  Parameter estimate(s):
    ##    lambda 
    ## 0.6655688 
    ## 
    ## Coefficients:
    ##                        Value  Std.Error   t-value p-value
    ## (Intercept)        0.2747425 0.29251762  0.939234  0.3479
    ## EFN                0.4743513 0.10338539  4.588186  0.0000
    ## Domatia            0.0217347 0.15247202  0.142549  0.8867
    ## Seed_Dispersal     0.2950486 0.07949651  3.711467  0.0002
    ## total.area.native  0.0000000 0.00000000  5.044367  0.0000
    ## abs_lat_native    -0.0141765 0.00259749 -5.457770  0.0000
    ## 
    ##  Correlation: 
    ##                   (Intr) EFN    Domati Sd_Dsp ttl.r.
    ## EFN               -0.011                            
    ## Domatia           -0.024 -0.138                     
    ## Seed_Dispersal     0.004 -0.112  0.022              
    ## total.area.native -0.053 -0.193  0.021 -0.286       
    ## abs_lat_native    -0.135  0.031  0.055 -0.005 -0.037
    ## 
    ## Standardized residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -1.28889508 -0.24969506 -0.03911948  0.20280009  4.05796633 
    ## 
    ## Residual standard error: 1.026786 
    ## Degrees of freedom: 769 total; 763 residual

``` r
## Diagnostic plots
plot(a3$residuals, a3$fitted)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-5.png)<!-- -->

``` r
qqnorm(a3$residuals)
qqline(a3$residuals)
```

![](README_files/figure-gfm/efns/domatia%20and%20ants%20models-6.png)<!-- -->
