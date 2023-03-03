---
title: "Generalized mutualism promotes range expansion in both ant and plant partners"
author: "Pooja Nathan and Megan Frederickson"
date: "`r Sys.Date()`"
output: github_document
editor_options: 
  chunk_output_type: console
---

# How does mutualism affect ant and plant range sizes?     
This R Markdown document describes the dataset and code for:

Nathan P, Economo EP, Guenard B, Simonsen A, Frederickson ME. Generalized mutualisms promote range expansion in both plant and ant partners. In prep. 

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

First, let's load the packages we'll use.

```{r packages, warning=FALSE, message=FALSE}
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
library(png)
library(patchwork)
library(grid)
```

## Legume dataset
We obtained introduced and native range size data for legumes from Simonsen et al. (2017).

Mutualistic trait data came from Weber et al. (2015) for extrafloral nectaries (EFNs), Chomicki & Renner (2015) for domatia, Simonsen et al. (2017) for nodulation, and Soudzilovskaia et al. (2020) for mycorrhizae. 

The next several chunks of code are slow to run, so they are not run here, but they are included for reproducibility.

```{r Read in mutualism and legume range datasets, warning=FALSE, message=FALSE, eval=FALSE}
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
We resolved misspelled names using the gnr_resolve function and checked for synonyms using the synonyms function in the taxize package. 

First, let's check the taxonomy of the names of domatia-bearing plants. 

```{r Check taxonomic names for domatia dataset, include=FALSE, warning=FALSE, message=FALSE, eval=FALSE}
#Not run
#Resolve names
domatia_resolve <- as.data.frame(gnr_resolve(domatia[, "Phy"], best_match_only=TRUE)) 
domatia <- merge(domatia, domatia_resolve, by.x = "Phy", by.y = "user_supplied_name", all.x=TRUE) #Merge resolved and original names
domatia$matched_name <- paste0(word(domatia$matched_name, 1, 1), " ", tolower(word(domatia$matched_name, 2, 2))) #Make lower case and reduce to just genus species (in other words, get rid of names of taxonomic authorities)

#Find synonyms for resolved names
domatia_synonyms_pow <- synonyms(domatia$matched_name, db="pow") #Get synonyms
domatia_syn <- as.data.frame(unlist(lapply(domatia_synonyms_pow, nth, 2))) #Make into dataframe
colnames(domatia_syn) <- "synonym" #Fix column name
domatia_syn$matched_name <- gsub('[[:digit:]]+', '', row.names(domatia_syn)) #Fix rownames

#Determine if original names and synonyms are in EFN, mycorrhizae, and range datasets
domatia_syn$synEFNY <- domatia_syn$synonym %in% EFN$Phy 
domatia_syn$synrangeY <- domatia_syn$synonym %in% range$Phy
domatia_syn$synmycoY <- domatia_syn$synonym %in% mycorrhizae$Phy
write.csv(domatia_syn, "domatia_synonyms.csv") #Save synonyms
domatia_syn <- subset(domatia_syn, domatia_syn$synEFNY | domatia_syn$synrangeY | domatia_syn$synmycoY) #Subset if synonyms match other datasets
domatia_syn$Domatia <- 1 #Add trait
colnames(domatia_syn)[[1]] <- "Phy" #Fix column name
rownames(domatia_syn) <- c() #Remove row names

#Add synonyms to domatia dataset
domatia <- rbind(domatia[, c("Domatia", "Phy", "matched_name")], domatia_syn[, c("Domatia", "Phy","matched_name")]) #Merge
domatia$Phy <- trimws(domatia$Phy) #Trim white space from taxonomic names

#Check if each domatia name is in the EFN, mycorrhizae, and range datasets
domatia$RangeY <- domatia$Phy %in% range$Phy
domatia$EFNY <- domatia$Phy %in% EFN$Phy
#domatia$mycoY <- domatia$Phy %in% mycorrhizae$Phy #No domatia-bearing plants in mycorrhizal dataset

write.csv(domatia, file="domatia_resolved.csv", row.names = FALSE)
```

Next, let's check the taxonomy of EFN-bearing plants.

```{r Check taxonomic names for EFN dataset, including=FALSE, warning=FALSE, message=FALSE, eval=FALSE}
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

```{r Merge EFN and domatia, warning=FALSE, message=FALSE}
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

Next, we'll check the taxonomy of the legumes for which we have data on mycorrhizae. 

```{r Check taxonomic names for mycorrhizae dataset, warning=FALSE, message=FALSE, eval=FALSE}
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

```{r Merge EFN, domatia, and mycorrhizae data, warning=FALSE, message=FALSE}
#Read in taxonomically resolved mycorrhizal dataset
mycorrhizae <- read.csv("mycorrhizae_resolved.csv")

#A little further cleaning of mycorrhizal dataset, to remove taxa that appear twice
mycorrhizae <- mycorrhizae %>% group_by(Phy, RangeY, domY, efnY) %>% dplyr::summarize(n.records=sum(n), sum.AM = sum(sum.AM), sum.EM=sum(sum.EM))
mycorrhizae$AM <- ifelse(mycorrhizae$sum.AM > 0, "Y", "N")
mycorrhizae$EM <- ifelse(mycorrhizae$sum.EM > 0, "Y", "N")

#Merge 
EFN_dom_myco <- merge(mycorrhizae, EFN_dom, by="Phy", all=TRUE)
```

Finally, we'll check the taxonomy of the range dataset.

```{r Check taxonomic names for range dataset, message=FALSE, warning=FALSE, eval=FALSE}
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

```{r Merge full legume datasets, warning=FALSE, message=FALSE}
range <- read.csv("range_resolved.csv")

legume_range_df <- merge(range, EFN_dom_myco, all.x="TRUE", all.y="FALSE", by= "Phy") #Put all the data in one dataframe
legume_range_df$EFN <- ifelse(is.na(legume_range_df$EFN), 0, legume_range_df$EFN) #Add zeros for NAs in EFN trait
legume_range_df$Domatia <- ifelse(is.na(legume_range_df$Domatia), 0, legume_range_df$Domatia) #Add zeros for NAs in domatia trait
write.csv(legume_range_df, file="legume_range_traits.csv", row.names = FALSE)
```

### Summarize legume dataset
We are finished cleaning the dataset, and now simply have to summarize, visualize, and model.

First, let's summarize how many species we have in each category. How many legumes with vs. without EFNs do we have range size data for, and how many introduced ranges have they been introduced to, on average? 

```{r Summarize legume EFN data in a table, warning=FALSE, message=FALSE}
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

#Assign to a shorter dataframe name
df <- legume_range_df

#Summarize the number of legume taxa with and without EFNs
summary.efn <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.efn)
```

How many legumes with vs. without domatia do we have range size data for? 

```{r Summarize legume domatia data in a table, warning=FALSE, message=FALSE}
summary.dom <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.dom)
```

How many legumes that do vs. do not form nodules do we have range size data for? 

```{r Summarize legume nodulation data in a table, warning=FALSE, message=FALSE}
summary.fix <- ungroup(subset(df, !is.na(num_introduced)) %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.fix)
```

How many legumes do vs. do not associate with mycorrhizae do we have range size data for?

```{r Summarize legume mycorrhizal data in a table, warning=FALSE, message=FALSE}
summary.myco <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(myco) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.myco)

summary.myco2 <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco)) %>% group_by(AM,EM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n), mean_area_introduced = mean(total_area_introduced, na.rm=TRUE), sd_area_introduced = sd(total_area_introduced, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n)))
kable(summary.myco2)
```

### Make figures
Let's make the figures showing how mutualistic trait affect legume native and introduced range sizes. 

```{r Legume plots for number of introduced ranges, warning=FALSE, message=FALSE, fig.width=12, fig.height=8, dpi=300}

#Set sizes and widths for all figures
pt_size <- 3
y_limits <- c(-0.5, 15)
er_width <- 0.1
y_text <- -0.25
y_inset_limits <- c(0,1)

df$introducedY <- ifelse(df$num_introduced > 0, 1, 0) #Create binary variable for whether or not legume is introduced

#EFN figure
summary.efn <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
#kable(summary.efn)

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
#kable(summary.dom)

p_dom <- ggplot(data=summary.dom, aes(x=Domatia, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "*")  

summary.dom2 <- ungroup(df %>% group_by(Domatia, introducedY) %>% dplyr::summarize(n=n()))
summary.dom2.wide <- spread(summary.dom2, key = introducedY, value=n)
colnames(summary.dom2.wide) <- c("Domatia","Not_introduced",  "Introduced")
summary.dom2.wide$total <- summary.dom2.wide$Not_introduced+summary.dom2.wide$Introduced
summary.dom2.wide$prop.introduced <- summary.dom2.wide$Introduced/summary.dom2.wide$total
prop.dom <- paste0(summary.dom2.wide$Introduced, "/", summary.dom2.wide$total)

inset_p_dom <- ggplot(data=summary.dom2.wide, aes(x=Domatia, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.dom), color="white")

#Nodules figure
summary.fix <- ungroup(subset(df, !is.na(num_introduced) & num_introduced > 0) %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
#kable(summary.fix)

p_fix <- ggplot(data=summary.fix, aes(x=fixer, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 13, label = "ns") 

summary.fix2 <- ungroup(df %>% group_by(fixer, introducedY) %>% dplyr::summarize(n=n()))
#kable(summary.fix2)
summary.fix2.wide <- spread(summary.fix2, key = introducedY, value=n)
colnames(summary.fix2.wide) <- c("fixer", "Not_introduced",  "Introduced")
summary.fix2.wide$total <- summary.fix2.wide$Not_introduced+summary.fix2.wide$Introduced
summary.fix2.wide$prop.introduced <- summary.fix2.wide$Introduced/summary.fix2.wide$total
prop.fix <- paste0(summary.fix2.wide$Introduced, "/", summary.fix2.wide$total)

inset_p_fix <- ggplot(data=summary.fix2.wide, aes(x=fixer, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Nodules")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "ns")+geom_text(aes(x=fixer, y=0.05, label=prop.fix), color="white")

#Mycorrhizae figure
summary.AM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(AM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
#kable(summary.AM)

p_AM <- ggplot(data=summary.AM, aes(x=AM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns")

summary.EM <- ungroup(subset(df, !is.na(num_introduced) & !is.na(myco) & num_introduced > 0) %>% group_by(EM) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n)))
kable(summary.EM)

p_EM <- ggplot(data=summary.EM, aes(x=EM, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EM, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+geom_line(aes(group=1), linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("EM")+geom_text(aes(x=EM, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 1.5, y = 13, label = "ns") 

summary.AM2 <- ungroup(subset(df, !is.na(myco)) %>% group_by(AM, introducedY) %>% dplyr::summarize(n=n()))
summary.AM2.wide <- spread(summary.AM2, key = introducedY, value=n)
colnames(summary.AM2.wide) <- c("AM", "Not_introduced",  "Introduced")
summary.AM2.wide$total <- summary.AM2.wide$Not_introduced+summary.AM2.wide$Introduced
summary.AM2.wide$prop.introduced <- summary.AM2.wide$Introduced/summary.AM2.wide$total
prop.AM <- paste0(summary.AM2.wide$Introduced, "/", summary.AM2.wide$total)

inset_p_AM <- ggplot(data=summary.AM2.wide, aes(x=AM, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("AM")+scale_y_continuous(limits=y_inset_limits)+annotate("text",x = 1.5, y = 0.55,label = "ns")+geom_text(aes(x=AM, y=0.05, label=prop.AM), color="white")

summary.EM2 <- ungroup(subset(df, !is.na(myco)) %>% group_by(EM, introducedY) %>% dplyr::summarize(n=n()))
summary.EM2.wide <- spread(summary.EM2, key = introducedY, value=n)
colnames(summary.EM2.wide) <- c("EM", "Not_introduced",  "Introduced")
summary.EM2.wide$total <- summary.EM2.wide$Not_introduced+summary.EM2.wide$Introduced
summary.EM2.wide$prop.introduced <- summary.EM2.wide$Introduced/summary.EM2.wide$total
prop.EM <- paste0(summary.EM2.wide$Introduced, "/", summary.EM2.wide$total)

inset_p_EM <- ggplot(data=summary.EM2.wide, aes(x=EM, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("EM")+scale_y_continuous(limits=y_inset_limits)+annotate("text",x = 1.5, y = 0.55,label = "ns")+geom_text(aes(x=EM, y=0.05, label=prop.EM), color="white")

fig1 <- plot_grid(inset_p_EFN, inset_p_dom, inset_p_fix, inset_p_AM, inset_p_EM, p_EFN, p_dom, p_fix, p_AM, p_EM, labels=c("AUTO"), nrow=2, ncol=5, align = "v", axis = "bt")
fig1
save_plot("Figure1.pdf", fig1, base_height=8, base_width=14)
```

We can also plot the effect of the same mutualistic traits on the native range size of legumes. 

```{r Legume plots for native area, warning=FALSE, message=FALSE, fig.width=8, fig.height=6, dpi=300}
pt_size <- 3
y_limits <- c(-500000, 13e+12)
er_width <- 0.1
y_text <- 0

#For EFNs
summary.efn.dom.native <- ungroup(df %>% group_by(EFN, Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_EFN_dom_native <- ggplot(data=summary.efn.dom.native, aes(x=EFN, y=mean_native, color=Domatia))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native, color=Domatia), width=er_width)+ geom_line(aes(group=Domatia), linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("EFNs")+geom_text(aes(x=c(0.8,1.2,1.8,2.2), y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+scale_color_grey(labels=c("No", "Yes"))+theme(legend.position = c(0.1, 0.8))+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")+annotate("text", x = 1.5, y = 3.9e+12, label = "**")

#For domatia
summary.dom.native <- ungroup(df %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_dom_native <- ggplot(data=summary.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

#For nodules
summary.fix.native <- ungroup(df %>% group_by(fixer) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))

p_fix_native <- ggplot(data=summary.fix.native, aes(x=fixer, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=fixer, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nodules")+geom_text(aes(x=fixer, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

#For EM fungi
summary.AM.native <- ungroup(df %>% group_by(AM) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.AM.native)

p_AM_native <- ggplot(data=subset(summary.AM.native, !is.na(AM)), aes(x=AM, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=AM, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("AM")+geom_text(aes(x=AM, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

#For EM fungi
summary.EM.native <- ungroup(df %>% group_by(EM) %>% dplyr::summarize(n=n(), mean_native = mean(total_area_native, na.rm=TRUE), sd_native = sd(total_area_native, na.rm=TRUE), se_native = sd_native/sqrt(n)))
kable(summary.EM.native)

p_EM_native <- ggplot(data=subset(summary.EM.native, !is.na(EM)), aes(x=EM, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=EM, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("EM")+geom_text(aes(x=EM, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = 1.5, y = 7.1e+12, label = "ns")

fig2 <- plot_grid(p_EFN_dom_native, p_fix_native, p_AM_native, p_EM_native, nrow=2, labels="AUTO", rel_widths = c(1, 1, 1, 1))
fig2
save_plot("Figure2.pdf", fig2, base_height = 8, base_width = 8)
```

### Statistical models

#### Mixed models
One approach is to use mixed models with legume tribe as a random effect to account for the non-independence of species in a tribe. The response variable (number of introduced ranges) is very non-normal so we fit two models: a binomial GLMM modelling whether or not a legume species is introduced, and linear mixed models of the number of introduced ranges for introduced species only. 

```{r Legume mixed models, warning=FALSE, message=FALSE}
#Fit binomial model for whether or not a legume species has been introduced
binomial1 <- glmer(introducedY~EFN+Domatia+fixer+uses_num_uses + +scale(total_area_native)+woody+annual+scale(abs_lat_native)+(1|tribe_ncbi), data=df, family="binomial",glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(binomial1) 
Anova(binomial1, type=3)
#plot(binomial1)

#Fit binomial model without random effect
binomial2 <- glm(introducedY~EFN+Domatia+fixer+uses_num_uses + +scale(total_area_native)+woody+annual+scale(abs_lat_native), data=df, family="binomial")
summary(binomial2) 
Anova(binomial2, type=3)
#plot(binomial2)

#Compare models
anova(binomial1, binomial2) #Random effect improves model fit to data

#Fit a linear mixed model for how many new ranges an introduced legume has established in, and how much total area they cover
legume_range_df_introducedY <- subset(df, num_introduced >0) #Filter to species with 1+ introduced ranges

#Number of introduced ranges
lmer1 <- lmer(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+woody+uses_num_uses+(1|tribe_ncbi), data=legume_range_df_introducedY)
summary(lmer1)
Anova(lmer1, type=3)
#plot(lmer1)

#Removed the interaction when nonsig 

#Number of introduced ranges without random effect of tribe
lm2 <- lm(log(num_introduced)~EFN+Domatia+fixer+scale(abs_lat_native)+scale(total_area_native)+annual+woody+uses_num_uses, data=legume_range_df_introducedY)
summary(lm2)
Anova(lm2, type=3)
#plot(lm2)

#Compare models
anova(lmer1, lm2) #Random effect improves model fit to data 

#Native range size
lmer3 <- lmer(log(total_area_native/1e+6)~EFN*Domatia+fixer+scale(abs_lat_native)+annual+woody+uses_num_uses+(1|tribe_ncbi), data=df)
summary(lmer3)
Anova(lmer3, type=3)


#Mycorrhizae
#Successful introduction?
binomial3 <- glmer(introducedY~AM+EM+scale(total_area_native)+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(df, !is.na(myco)), family="binomial", glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=100000)))
summary(binomial3) 
Anova(binomial3, type=3)
#plot(binomial3)

#Number of introduced ranges
lmer4 <- lmer(log(num_introduced)~EM+AM+scale(total_area_native)+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(legume_range_df_introducedY, !is.na(myco)))
summary(lmer4)
Anova(lmer4, type=3)
#plot(lmer4)

#Native range size
lmer5 <- lmer(log(total_area_native/1e+6)~AM+EM+annual+woody+scale(abs_lat_native)+uses_num_uses+(1|tribe_ncbi), data=subset(df, !is.na(myco)))
summary(lmer5)
Anova(lmer5, type=3)
#plot(lmer5)
```

#### Preparing dataset for pgls analysis


```{r PGLS, warning=FALSE, message=FALSE}
#Creating a column with genus and species name in the appropriate format
df$Phy2 <- paste0(as.character(word(legume_range_df$Phy, 1, 1)), "_", as.character(word(legume_range_df$Phy, 2, 2))) #sfsg

zanne <- read.tree("Vascular_Plants_rooted.dated.tre") #reading in Zanne et al. 2014 plant phylogeny
phyint <- intersect(zanne$tip.label, df$Phy2)  
phydiff <- setdiff(zanne$tip.label, df$Phy2)
pruned.tree.pgls <- drop.tip(zanne, phydiff) #dropping tips not in the dataset

range_pgls <- df[df$Phy2 %in% phyint, ]
#colnames(range_pgls)
#which(colSums(is.na(range_pgls))>0) #Check which columns have NAs
range_pgls <-range_pgls[,-c(23:30, 32, 34)] #Remove some unneeded columns
range_pgls <- range_pgls[complete.cases(range_pgls), ] #removing NA elements
```

#### PGLS models for EFN, domatia, and rhizobia

We ran PGLS models for legumes to account for phylogenetic relationships when measuring the effects of mutualisms on introduction success and introduced range size. The problem with PGLS models is that we lose a lot of trait data for species not in the phylogeny. We either need a better phylogeny or a non-phylogenetic model. 

```{r Legume pgls models for EFN domatia and rhizobia, warning=FALSE, message=FALSE}
#PGLS of number of introduced ranges as response variable with EFNs, domatia, and nodules, 
#interaction and other covariates
#Interactions that were not significant were removed

pgls1 <- gls(log(num_introduced + 1) ~ EFN+fixer+Domatia+total_area_native+abs_lat_native+woody+uses_num_uses+annual, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method="ML", data = range_pgls) 
summary(pgls1)


#Repeating for native area
#PGLS with both EFN and domatia presence and fixer, interaction and covariates
pgls2 <- gls(log((total_area_native/1e+6) + 1) ~ EFN*Domatia + fixer+ abs_lat_native+ annual + woody + uses_num_uses, correlation = corPagel(1, phy = pruned.tree.pgls, form = ~ Phy2), method = "ML", data = range_pgls)
summary(pgls2)

```

#### PGLS models for mycorrhizae

```{r Legume pgls models for mycorrhizae, warning=FALSE, message=FALSE }
range_myco <- subset(df, !is.na(myco))
range_myco$Phy2 <- as.character(range_myco$Phy2)
phyint1 <- intersect(zanne$tip.label, range_myco$Phy2)  
phydiff1 <- setdiff(zanne$tip.label, range_myco$Phy2)
pruned.myco.pgls <- drop.tip(zanne, phydiff1) #dropping tips not in the dataset

range_myco_pgls <- range_myco[range_myco$Phy2 %in% phyint1, ]
#which(colSums(is.na(range_myco_pgls))>0) #Check which columns have NAs
range_myco_pgls <-range_myco_pgls[,-c(32)] #Remove some unneeded columns
range_myco_pgls <- range_myco_pgls[complete.cases(range_myco_pgls), ] #removing NA elements

#PGLS of number of introduced ranges as response variable mycorrhizae and covariates as predictors
pgls3 <- gls(log(num_introduced + 1) ~ AM+EM + total_area_native + abs_lat_native + annual + uses_num_uses +woody, correlation = corPagel(1, phy = pruned.myco.pgls, form = ~ Phy2), method = "ML", data = range_myco_pgls) 
summary(pgls3)

#Repeating for native area
pgls4 <- gls((total_area_native)^(1/4) ~ AM+EM+ abs_lat_native+ annual + woody + uses_num_uses, correlation = corPagel(1, phy = pruned.myco.pgls, form = ~ Phy2), method = "ML", data = range_myco_pgls)
summary(pgls4)
```

#### Multiple mutualisms

We can also look at the effects of multiple mutualisms on introduction success and no. of introduced ranges in legumes.

```{r Number of mutualisms, warning=FALSE, message=FALSE, fig.width=6, fig.height=4, dpi=300}
#Creating a numeric vector for AM and EM states
df$AM.num <- ifelse(df$AM == "Y", 1, ifelse(df$AM == "N", 0, NA))
df$EM.num <- ifelse(df$EM == "Y", 1, ifelse(df$EM == "N", 0, NA))

#Creating a column containing the total number of mutualisms that each species participates in
df$num_mutualisms <- as.numeric(as.character(df$fixer))+as.numeric(as.character(df$EFN))+as.numeric(as.character(df$AM.num))+as.numeric(as.character(df$Domatia))+as.numeric(as.character(df$EM.num))

#Making figures for proportion of species introduced
summary.mnum2 <- ungroup(subset(df, !is.na(num_mutualisms)) %>% group_by(num_mutualisms, introducedY) %>% dplyr::summarize(n=n()))
summary.mnum2.wide <- spread(summary.mnum2, key = introducedY, value=n)
colnames(summary.mnum2.wide) <- c("num_mutualisms","Not_introduced",  "Introduced")
summary.mnum2.wide$total <- summary.mnum2.wide$Not_introduced+summary.mnum2.wide$Introduced
summary.mnum2.wide$prop.introduced <- summary.mnum2.wide$Introduced/(summary.mnum2.wide$total)
prop.mnum <- paste0(summary.mnum2.wide$Introduced, "/", summary.mnum2.wide$total)

inset_p_mnum <- ggplot(data=subset(summary.mnum2.wide, !is.na(num_mutualisms)), aes(x=num_mutualisms, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+xlab("Mutualisms (no.)")+ylab("Introduced (prop.)")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 2, y = 0.55, label = "*")+geom_text(aes(x=num_mutualisms, y=0.05, label=prop.mnum), color="white")

#Making figures for number of introduced ranges
summary.mnum <- subset(df, !is.na(num_mutualisms) & introducedY == 1) %>% group_by(num_mutualisms) %>% dplyr::summarize(n=n(), mean_num_introduced = mean(num_introduced, na.rm=TRUE), sd_num_introduced = sd(num_introduced, na.rm=TRUE), se_num_introduced = sd_num_introduced/sqrt(n))
summary.mnum$num_mutualisms <- as.factor(summary.mnum$num_mutualisms)

p_num_mutualisms <- ggplot(data=summary.mnum, aes(x=num_mutualisms, y=mean_num_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num_mutualisms, ymin=mean_num_introduced-se_num_introduced, ymax=mean_num_introduced+se_num_introduced), width=er_width)+theme_cowplot()+ylab("Introduced ranges (no.)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num_mutualisms, y= 3, label=n))+annotate("text", x=3, y=15.5, label="p = 0.07")
p_num_mutualisms

#Linear models to see if the proportion introduced and the number of introduced ranges vary with number of mutualisms
lmer6 <-  lmer(log(num_introduced)~num_mutualisms+annual+scale(abs_lat_native)+woody+uses_num_uses+scale(total_area_native)+(1|tribe_ncbi), data=subset(df, !is.na(num_mutualisms) & introducedY == 1))
summary(lmer6)
Anova(lmer6, type=3)
#plot(lmer6)

binomial4 <-  glmer(introducedY~num_mutualisms+annual+woody+scale(abs_lat_native)+uses_num_uses+scale(total_area_native)+(1|tribe_ncbi), data=subset(df, !is.na(num_mutualisms)), family="binomial")
summary(binomial4)
Anova(binomial4, type=3)
#plot(binomial4)

fig3 <- plot_grid(inset_p_mnum, p_num_mutualisms, ncol=1, labels="AUTO")
fig3
save_plot("Figure3.pdf", fig3, base_height = 8, base_width=5)

sum.nmum <- subset(df, !is.na(num_mutualisms)) %>% group_by(AM, EM, Domatia, EFN, fixer, num_mutualisms) %>% dplyr::summarize(n=n())
kable(sum.nmum)
```

#### Pagel correlations between traits
Pagel's correlation is used to understand if two binary traits show correlated evolution on a phylogeny.
Here, we calculated Pagel correlation for each pair of mutualistic traits in legumes.


```{r Legume Pagel, warning=FALSE, message=FALSE, eval=FALSE}
#Pagel correlation between traits
range_pagel <- df[df$Phy2 %in% phyint, ]
range_pagel <-range_pagel[,-c(23:28, 32, 34)] #Remove some unneeded columns

range_pagel <- range_pagel[!duplicated(range_pagel$Phy2), ] #Removing duplicated data
range_pagel <- range_pagel[complete.cases(range_pagel), ] 

names <- pruned.tree.pgls$tip.label 
names <- names[names %in% range_pagel$Phy2]

traits <- data.frame(range_pagel$Phy2, range_pagel$AM, range_pagel$EM, range_pagel$EFN, range_pagel$Domatia, range_pagel$fixer)
colnames(traits) <- c("Phy", "AM", "EM", "EFN", "Domatia", "Fixer")
traits <- arrange(traits, names)

#[match(as.character(pruned.tree.pgls$tip.label), traits$Phy),]

rownames(traits) <- traits[,1]
traits[,1] <- NULL
head(traits)

AM<-setNames(traits$AM,rownames(traits))
EM<-setNames(traits$EM,rownames(traits))
EFN<-setNames(traits$EFN,rownames(traits))
Domatia<-setNames(traits$Domatia,rownames(traits))
fixer<-setNames(traits$Fixer,rownames(traits))

#Subsetting phylogeny
phyint2 <- intersect(pruned.tree.pgls$tip.label, range_pagel$Phy2)  
phydiff2 <- setdiff(pruned.tree.pgls$tip.label, range_pagel$Phy2)
pruned.pagel <- drop.tip(pruned.tree.pgls, phydiff2) #dropping tips not in the dataset

fit.amem <- fitPagel(pruned.pagel, AM, EM)
amem <- plot(fit.amem,lwd.by.rate=TRUE)

fit.amefn <- fitPagel(pruned.pagel, AM, EFN) #nonsig

fit.amdom <- fitPagel(pruned.pagel, AM, Domatia) #signif.
amdom <- plot(fit.amdom,lwd.by.rate=TRUE)

fit.amfix <- fitPagel(pruned.pagel, AM, fixer) #pval 0.009
amfix <- plot(fit.amfix,lwd.by.rate=TRUE)

fit.emefn <- fitPagel(pruned.pagel, EM, EFN) #nonsig

fit.emdom <- fitPagel(pruned.pagel, EM, Domatia)
emdom <- plot(fit.emdom,lwd.by.rate=TRUE)

fit.emfix <- fitPagel(pruned.pagel, EM, fixer) #pval 0.0004
emfix <- plot(fit.emfix, lwd.by.rate=TRUE)

fit.efndom <- fitPagel(pruned.pagel, EFN, Domatia)
efndom <- plot(fit.efndom, lwd.by.rate=TRUE)

fit.efnfix <- fitPagel(pruned.pagel, EFN, fixer) #nonsig

fit.domfix <- fitPagel(pruned.pagel, Domatia, fixer) #marginally sig
domfix <- plot(fit.domfix, lwd.by.rate=TRUE)

```

## Ant dataset

We obtained introduced and native range size data for ants from Benoit Guenard and Evan Economo, and trait data on which ant species visit EFNs, disperse seeds, and nest in domatia from Kaur et al. (2019).

```{r efns and ants, warning=FALSE, message=FALSE}
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

We will use the ant tribe as a random effect in models, and so need to download it for each ant genus in the dataset.

```{r Get ant tribe, eval=FALSE, warning=FALSE, message=FALSE}
#Not run
#Code included for reproducibility
ant <- data.frame(genus = unique(word(area$Phy, 1, 1, sep="_"))) #creating the vector of unique genera
ant$tribe <- NA #create an empty vector ant tribe 
#Loop through all ant genera
for (i in 1:length(ant$genus)){
  ant[i, "tribe"] <- tax_name(sci= ant$genus[i], get = 'tribe', db = 'itis')$tribe
}
write.csv(ant, "ant_tribe.csv") 
#some tribes added manually from AntWiki because data was not available

ant <- read.csv("ant_tribe.csv")
#Merge ant tribe with ant range dataset
area$genus <- word(area$Phy, 1, 1, sep="_")
area <- merge(area, ant, by="genus", all.x=TRUE)
write.csv(area, "ant_areas.csv")
```

### Summarize ant dataset

```{r}
area <- read.csv("ant_areas.csv")
area$tribe <- as.factor(area$tribe)

sum(complete.cases(subset(area, !is.na(introducedY))[, c("EFN", "Domatia", "Seed_Dispersal")])) #Total taxa with at least some trait data 
kable(subset(area, !is.na(introducedY)) %>% group_by(EFN) %>% dplyr::summarize(n=n()))
kable(subset(area, !is.na(introducedY)) %>% group_by(Domatia) %>% dplyr::summarize(n=n()))
kable(subset(area, !is.na(introducedY)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n()))
```
### Make figures
Let's make the figures showing how mutualistic trait affect ant native and introduced range sizes. 

```{r Ant plots for number of introduced ranges, warning=FALSE, message=FALSE}
pt_size <- 3
y_limits <- c(3, 36)
er_width <- 0.1
y_text <- 4


#For EFN visitation
summary.ant.efn <- subset(area, introducedY == 1 & !is.na(EFN)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antEFN <- ggplot(data=summary.ant.efn, aes(x=EFN, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Visits EFNs")+scale_x_discrete(labels=c("No", "Yes"))+geom_text(aes(x=EFN, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+annotate("text", x = 0.5, y = 20, label = "**")
#Adjust labels

summary.ant.efn2 <- ungroup(subset(area, !is.na(EFN)) %>% group_by(EFN, introducedY) %>% dplyr::summarize(n=n()))
summary.ant.efn2.wide <- spread(summary.ant.efn2, key = introducedY, value=n)
colnames(summary.ant.efn2.wide) <- c("EFN","Not_introduced",  "Introduced")
summary.ant.efn2.wide$total <- summary.ant.efn2.wide$Not_introduced+summary.ant.efn2.wide$Introduced
summary.ant.efn2.wide$prop.introduced <- summary.ant.efn2.wide$Introduced/(summary.ant.efn2.wide$total)
prop.ant.efn <- paste0(summary.ant.efn2.wide$Introduced, "/", summary.ant.efn2.wide$total)

inset_p_antEFN <- ggplot(data=summary.ant.efn2.wide, aes(x=EFN, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 0.5, y = 0.55, label = "***")+geom_text(aes(x=EFN, y=0.05, label=prop.ant.efn), color="white")

#For residing in domatia
summary.ant.dom <- subset(area, introducedY == 1 & !is.na(Domatia)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_antdom <- ggplot(data=summary.ant.dom, aes(x=Domatia, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Nests in domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 0.5, y = 20, label = "ns")

summary.ant.dom2 <- ungroup(subset(area, !is.na(Domatia)) %>% group_by(Domatia, introducedY) %>% dplyr::summarize(n=n()))
summary.ant.dom2.wide <- spread(summary.ant.dom2, key = introducedY, value=n)
colnames(summary.ant.dom2.wide) <- c("Domatia","Not_introduced",  "Introduced")
summary.ant.dom2.wide$total <- summary.ant.dom2.wide$Not_introduced+summary.ant.dom2.wide$Introduced
summary.ant.dom2.wide$prop.introduced <- summary.ant.dom2.wide$Introduced/(summary.ant.dom2.wide$total)
prop.ant.dom <- paste0(summary.ant.dom2.wide$Introduced, "/", summary.ant.dom2.wide$total)

inset_p_antdom <- ggplot(data=summary.ant.dom2.wide, aes(x=Domatia, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 0.5, y = 0.55, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.ant.dom), color="white")

#For myrmecochory (seed dispersal by ants)
summary.ant.elaiosome <- subset(area, introducedY == 1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n(), mean_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_introduced = sd_introduced/sqrt(n))

p_elaiosome <- ggplot(data=summary.ant.elaiosome, aes(x=Seed_Dispersal, y=mean_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=Seed_Dispersal, ymin=mean_introduced-se_introduced, ymax=mean_introduced+se_introduced), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Introduced ranges (no.)")+xlab("Disperses seeds")+geom_text(aes(x=Seed_Dispersal, y= y_text, label=n))+scale_y_continuous(limits=y_limits)+scale_x_discrete(labels=c("No", "Yes"))+annotate("text", x = 0.5, y = 20, label = "ns")

summary.ant.seed2 <- ungroup(subset(area, !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, introducedY) %>% dplyr::summarize(n=n()))
summary.ant.seed2.wide <- spread(summary.ant.seed2, key = introducedY, value=n)
colnames(summary.ant.seed2.wide) <- c("Seed_Dispersal","Not_introduced",  "Introduced")
summary.ant.seed2.wide$total <- summary.ant.seed2.wide$Not_introduced+summary.ant.seed2.wide$Introduced
summary.ant.seed2.wide$prop.introduced <- summary.ant.seed2.wide$Introduced/(summary.ant.seed2.wide$total)
prop.ant.seed <- paste0(summary.ant.seed2.wide$Introduced, "/", summary.ant.seed2.wide$total)

inset_p_antseed <- ggplot(data=summary.ant.seed2.wide, aes(x=Seed_Dispersal, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+scale_x_discrete(labels=c("No", "Yes"))+ylab("Introduced (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 0.5, y = 0.55, label = "***")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.ant.seed), color="white")

fig4 <- plot_grid(inset_p_antEFN, inset_p_antseed, inset_p_antdom, p_antEFN, p_elaiosome, p_antdom, nrow=2, labels="AUTO", align="v", axis="bt")
fig4
save_plot("Figure4.pdf", fig4, base_width=14, base_height =8)
```

#### Native range size
Making figures to show the relationship between native range and mutualistic interactions in ants.

```{r}
y_limits = c(3200000, 1.25e+07)
y_text = 3300000

#For visiting EFNs
summary.ant.efn.native <- subset(area, !is.na(EFN) & !is.na(total.area.native)) %>% group_by(EFN) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))


p_antEFN_native <- ggplot(data=summary.ant.efn.native, aes(x=EFN, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=EFN, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Visits EFNs")+geom_text(aes(x=EFN, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = .5, y = 8100000, label = "***")

#For resifding in domatia
summary.ant.dom.native <- subset(area, !is.na(Domatia)& !is.na(total.area.native)) %>% group_by(Domatia) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antdom_native <- ggplot(data=summary.ant.dom.native, aes(x=Domatia, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Domatia, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Nests in domatia")+geom_text(aes(x=Domatia, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = .5, y = 8100000, label = "ns")

#For seed dispersal
summary.ant.seed.native <- subset(area, !is.na(Seed_Dispersal)& !is.na(total.area.native)) %>% group_by(Seed_Dispersal) %>% dplyr::summarize(n=n(), mean_native = mean(total.area.native, na.rm=TRUE), sd_native = sd(total.area.native, na.rm=TRUE), se_native = sd_native/sqrt(n))

p_antseed_native <- ggplot(data=summary.ant.seed.native, aes(x=Seed_Dispersal, y=mean_native))+geom_point(size=pt_size)+geom_errorbar(aes(x=Seed_Dispersal, ymin=mean_native-se_native, ymax=mean_native+se_native), width=er_width)+ geom_line(aes(group=1),linetype="dashed")+theme_cowplot()+ylab("Native range (sq. km)")+xlab("Disperses seeds")+geom_text(aes(x=Seed_Dispersal, y= y_text, label=n))+scale_x_discrete(labels=c("No", "Yes"))+scale_y_continuous(limits=y_limits)+annotate("text", x = .5, y = 8100000, label = "***")

fig5 <- plot_grid(p_antEFN_native, p_antseed_native, p_antdom_native, nrow=1, labels="AUTO")
fig5
save_plot("Figure5.pdf", fig5, base_height = 4, base_width = 12)
save_plot("Figure5.png", fig5, base_height = 4, base_width = 12)
```

### Mixed models
We used the same approach for analyzing the ant dataset as we did for the legumes. The response variable (number of introduced ranges) is very non-normal so we fit two models: a binomial GLMM modelling whether or not a legume species is introduced, and linear mixed models of the number of introduced ranges for introduced species only. We also used Tribe as a random effect.

```{r ant glmms, warning=FALSE, message=FALSE}
#Checking the dependence of the proportion of species introduced and the number of non-contiguous ranges on mutualism

#Introduction success
binomial5 <- lme4::glmer(introducedY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=area, family="binomial")
summary(binomial5)
car::Anova(binomial5, type=3)
#plot(binomial5)

#Total non-contiguous ranges
lmer7 <- lmer(log(n_introduced_ranges)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=subset(area, introducedY == 1))
summary(lmer7)
Anova(lmer7, type=3)
#plot(lmer7)

#Native range area
lmer10 <- lmer(log(total.area.native)~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+(1|tribe), data=area)
summary(lmer10)
Anova(lmer10, type=3)
#plot(lmer10)
```


### Multiple mutualisms
We can also look at the effects of multiple mutualisms on introduction success and no. of introduced ranges in ants


```{r}
#Calculating the number of mutualisms that each ant species participates in
area$num.mm <- as.numeric(as.character(area$Seed_Dispersal))+as.numeric(as.character(area$Domatia))+as.numeric(as.character(area$EFN))

#Model to understand if introduction success varies with number of mutualisms
binomial10 <- glmer(introducedY ~ num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=area, family ="binomial")
summary(binomial10)
Anova(binomial10, type=3)
plot(binomial10)

#Model to understand if the number of introduced ranges vary with number of mutualisms
lmer11 <- lmer(log(n_introduced_ranges)~num.mm + scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=subset(area, introducedY == 1))
summary(lmer11)
Anova(lmer11, type=3)

#Making figures to understand the effects of the number of mutualisms on introduction success
#and number of introiduced ranges

#Introduction success
summary.ant.mm <- ungroup(subset(area, !is.na(introducedY)) %>% group_by(num.mm, introducedY) %>% dplyr::summarize(n=n()))
summary.ant.mm.wide <- spread(summary.ant.mm, key=introducedY, value=n)
colnames(summary.ant.mm.wide) <- c("Num_mutualisms","Not_introduced",  "Introduced")
summary.ant.mm.wide$total <- summary.ant.mm.wide$Not_introduced+summary.ant.mm.wide$Introduced
summary.ant.mm.wide$prop.introduced <- summary.ant.mm.wide$Introduced/(summary.ant.mm.wide$total)
prop.ant.mm <- paste0(summary.ant.mm.wide$Introduced, "/", summary.ant.mm.wide$total)
prop.ant.mm <- prop.ant.mm[1:4]

inset_p_antmm <- ggplot(data=subset(summary.ant.mm.wide, !is.na(Num_mutualisms)),aes(x=Num_mutualisms, y=prop.introduced))+geom_bar(stat="identity")+theme_cowplot()+ylab("Introduced (prop.)")+xlab("Mutualisms (no.)")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.55, label = "***")+geom_text(aes(x=Num_mutualisms, y=0.05, label=prop.ant.mm), color="white")

#Number of introduced ranges
summary.ant.mnum1 <- subset(area, !is.na(num.mm) & introducedY == 1) %>% group_by(num.mm) %>% dplyr::summarize(n=n(), mean_area_introduced = mean(n_introduced_ranges, na.rm=TRUE), sd_area_introduced = sd(n_introduced_ranges, na.rm=TRUE), se_area_introduced = sd_area_introduced/sqrt(n))
summary.ant.mnum1$num.mm <- as.factor(summary.ant.mnum1$num.mm)

p_num_mm <- ggplot(data=summary.ant.mnum1, aes(x=num.mm, y=mean_area_introduced))+geom_point(size=pt_size)+geom_errorbar(aes(x=num.mm, ymin=mean_area_introduced-se_area_introduced, ymax=mean_area_introduced+se_area_introduced), width=er_width)+theme_cowplot()+ylab("Introduced range (sq. km)")+geom_line(aes(group=1),linetype="dashed")+xlab("Mutualisms (no.)")+geom_text(aes(x=num.mm, y= -1, label=n))+annotate("text", x=1.5, y=25, label="**")

fig6 <- plot_grid(inset_p_antmm, p_num_mm, ncol=1, labels="AUTO")
fig6
save_plot("Figure6.pdf", fig6, base_height = 8, base_width=5)
```

### By introduction mode
We also have data on the introduction mode of ants: in their new ranges, ant species may have been naturalized (exotic), introduced indoors, or intercepted by customs. In this analysis, we trid to understand how participation in mutualisms affected the mode of introduction.

```{r Introduction mode}
#Making figures 

# Exotic
area$ExoticY <- as.factor(area$ExoticY)
exotic.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, ExoticY) %>% dplyr::summarize(n=n()))
exotic.efn.wide <- spread(exotic.efn, key = ExoticY, value=n)
colnames(exotic.efn.wide) <- c("EFN","Not_exotic",  "Exotic")
exotic.efn.wide$EFN <- c("No", "Yes")
exotic.efn.wide$total <- exotic.efn.wide$Not_exotic + exotic.efn.wide$Exotic
exotic.efn.wide$prop.exotic <- exotic.efn.wide$Exotic/exotic.efn.wide$total
prop.exotic.efn <- paste0(exotic.efn.wide$Exotic, "/", exotic.efn.wide$total)

p_exoticefn <- ggplot(data=exotic.efn.wide, aes(x=EFN, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=EFN, y=0.05, label=prop.exotic.efn), color="white")

exotic.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, ExoticY) %>% dplyr::summarize(n=n()))
exotic.dom.wide <- spread(exotic.dom, key = ExoticY, value=n)
colnames(exotic.dom.wide) <- c("Domatia","Not_exotic",  "Exotic")
exotic.dom.wide$Domatia <- c("No", "Yes")
exotic.dom.wide$total <- exotic.dom.wide$Not_exotic + exotic.dom.wide$Exotic
exotic.dom.wide$prop.exotic <- exotic.dom.wide$Exotic/exotic.dom.wide$total
prop.exotic.dom <- paste0(exotic.dom.wide$Exotic, "/", exotic.dom.wide$total)

p_exoticdom <- ggplot(data=exotic.dom.wide, aes(x=Domatia, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.exotic.dom), color="white")

exotic.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, ExoticY) %>% dplyr::summarize(n=n()))
exotic.seed.wide <- spread(exotic.seed, key = ExoticY, value=n)
colnames(exotic.seed.wide) <- c("Seed_Dispersal","Not_exotic",  "Exotic")
exotic.seed.wide$Seed_Dispersal <- c("No", "Yes")
exotic.seed.wide$total <- exotic.seed.wide$Not_exotic + exotic.seed.wide$Exotic
exotic.seed.wide$prop.exotic <- exotic.seed.wide$Exotic/exotic.seed.wide$total
prop.exotic.seed <- paste0(exotic.seed.wide$Exotic, "/", exotic.seed.wide$total)

p_exoticseed <- ggplot(data=exotic.seed.wide, aes(x=Seed_Dispersal, y=prop.exotic))+geom_bar(stat="identity")+theme_cowplot()+ylab("Naturalized (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.exotic.seed), color="white")

p_exotic <- plot_grid(p_exoticefn, p_exoticdom, p_exoticseed, nrow=1, labels=c("G", "H", "I"))

#Binomial model to understand how the introduction success of naturalized species depends of mutualism participation
binomial7 <- glmer(ExoticY~EFN+Domatia+Seed_Dispersal+scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=subset(area, introducedY == 1), family="binomial")
summary(binomial7)
Anova(binomial7, type=3)
#plot(binomial7)

#Indoors?
area$IndoorY <- as.factor(area$IndoorY)

indoor.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, IndoorY) %>% dplyr::summarize(n=n()))
indoor.efn.wide <- spread(indoor.efn, key = IndoorY, value=n)
colnames(indoor.efn.wide) <- c("EFN","Not_indoor",  "indoor")
indoor.efn.wide$EFN <- c("No", "Yes")
indoor.efn.wide$total <- indoor.efn.wide$Not_indoor + indoor.efn.wide$indoor
indoor.efn.wide$prop.indoor <- indoor.efn.wide$indoor/indoor.efn.wide$total
prop.indoor.efn <- paste0(indoor.efn.wide$indoor, "/", indoor.efn.wide$total)

p_indoorefn <- ggplot(data=indoor.efn.wide, aes(x=EFN, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "**")+geom_text(aes(x=EFN, y=0.05, label=prop.indoor.efn), color="white")

indoor.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, IndoorY) %>% dplyr::summarize(n=n()))
indoor.dom.wide <- spread(indoor.dom, key = IndoorY, value=n)
colnames(indoor.dom.wide) <- c("Domatia","Not_indoor",  "indoor")
indoor.dom.wide$Domatia <- c("No", "Yes")
indoor.dom.wide$total <- indoor.dom.wide$Not_indoor + indoor.dom.wide$indoor
indoor.dom.wide$prop.indoor <- indoor.dom.wide$indoor/indoor.dom.wide$total
prop.indoor.dom <- paste0(indoor.dom.wide$indoor, "/", indoor.dom.wide$total)

p_indoordom <- ggplot(data=indoor.dom.wide, aes(x=Domatia, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.indoor.dom), color="white")

indoor.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, IndoorY) %>% dplyr::summarize(n=n()))
indoor.seed.wide <- spread(indoor.seed, key = IndoorY, value=n)
colnames(indoor.seed.wide) <- c("Seed_Dispersal","Not_indoor",  "indoor")
indoor.seed.wide$Seed_Dispersal <- c("No", "Yes")
indoor.seed.wide$total <- indoor.seed.wide$Not_indoor + indoor.seed.wide$indoor
indoor.seed.wide$prop.indoor <- indoor.seed.wide$indoor/indoor.seed.wide$total
prop.indoor.seed <- paste0(indoor.seed.wide$indoor, "/", indoor.seed.wide$total)

p_indoorseed <- ggplot(data=indoor.seed.wide, aes(x=Seed_Dispersal, y=prop.indoor))+geom_bar(stat="identity")+theme_cowplot()+ylab("Indoors (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "**")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.indoor.seed), color="white")

p_indoor <- plot_grid(p_indoorefn, p_indoordom, p_indoorseed, nrow=1, labels=c("D","E", "F"))

#Binomial model to understand how the introduction success of indoor introduced species depends of mutualism participation

binomial8 <- glmer(IndoorY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=subset(area, introducedY == 1), family="binomial")
summary(binomial8)
Anova(binomial8, type=3)
#plot(binomial8)

#Intercepted
area$InterceptedY <- as.factor(area$InterceptedY)

intercepted.efn <- ungroup(subset(area, introducedY ==1 & !is.na(EFN)) %>% group_by(EFN, InterceptedY) %>% dplyr::summarize(n=n()))
intercepted.efn.wide <- spread(intercepted.efn, key = InterceptedY, value=n)
colnames(intercepted.efn.wide) <- c("EFN","Not_intercepted",  "intercepted")
intercepted.efn.wide$EFN <- c("No", "Yes")
intercepted.efn.wide$total <- intercepted.efn.wide$Not_intercepted + intercepted.efn.wide$intercepted
intercepted.efn.wide$prop.intercepted <- intercepted.efn.wide$intercepted/intercepted.efn.wide$total
prop.intercepted.efn <- paste0(intercepted.efn.wide$intercepted, "/", intercepted.efn.wide$total)

p_interceptedefn <- ggplot(data=intercepted.efn.wide, aes(x=EFN, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Visits EFNs")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "ns")+geom_text(aes(x=EFN, y=0.05, label=prop.intercepted.efn), color="white")

intercepted.dom <- ungroup(subset(area, introducedY ==1 & !is.na(Domatia)) %>% group_by(Domatia, InterceptedY) %>% dplyr::summarize(n=n()))
intercepted.dom.wide <- spread(intercepted.dom, key = InterceptedY, value=n)
colnames(intercepted.dom.wide) <- c("Domatia","Not_intercepted",  "intercepted")
intercepted.dom.wide$Not_intercepted <- ifelse(is.na(intercepted.dom.wide$Not_intercepted), 0, intercepted.dom.wide$Not_intercepted)
intercepted.dom.wide$Domatia <- c("No", "Yes")
intercepted.dom.wide$total <- intercepted.dom.wide$Not_intercepted + intercepted.dom.wide$intercepted
intercepted.dom.wide$prop.intercepted <- intercepted.dom.wide$intercepted/intercepted.dom.wide$total
prop.intercepted.dom <- paste0(intercepted.dom.wide$intercepted, "/", intercepted.dom.wide$total)

p_intercepteddom <- ggplot(data=intercepted.dom.wide, aes(x=Domatia, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Nests in domatia")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.45, y = 0.85, label = "ns")+geom_text(aes(x=Domatia, y=0.05, label=prop.intercepted.dom), color="white")

intercepted.seed <- ungroup(subset(area, introducedY ==1 & !is.na(Seed_Dispersal)) %>% group_by(Seed_Dispersal, InterceptedY) %>% dplyr::summarize(n=n()))
intercepted.seed.wide <- spread(intercepted.seed, key = InterceptedY, value=n)
colnames(intercepted.seed.wide) <- c("Seed_Dispersal","Not_intercepted",  "intercepted")
intercepted.seed.wide$Seed_Dispersal <- c("No", "Yes")
intercepted.seed.wide$total <- intercepted.seed.wide$Not_intercepted + intercepted.seed.wide$intercepted
intercepted.seed.wide$prop.intercepted <- intercepted.seed.wide$intercepted/intercepted.seed.wide$total
prop.intercepted.seed <- paste0(intercepted.seed.wide$intercepted, "/", intercepted.seed.wide$total)

p_interceptedseed <- ggplot(data=intercepted.seed.wide, aes(x=Seed_Dispersal, y=prop.intercepted))+geom_bar(stat="identity")+theme_cowplot()+ylab("Intercepted (prop.)")+xlab("Disperses seeds")+scale_y_continuous(limits=y_inset_limits)+annotate("text", x = 1.5, y = 0.85, label = "***")+geom_text(aes(x=Seed_Dispersal, y=0.05, label=prop.intercepted.seed), color="white")

p_intercepted <- plot_grid(p_interceptedefn, p_intercepteddom, p_interceptedseed, nrow=1, labels="AUTO")

#Binomial model to understand how the introduction success of species intercepted by customs depends of mutualism participation

binomial9 <- glmer(InterceptedY~EFN+Seed_Dispersal+Domatia+scale(abs_lat_native)+scale(total.area.native)+(1|tribe), data=subset(area, introducedY == 1), family="binomial")
summary(binomial9)
Anova(binomial9, type=3)
#plot(binomial9)

fig7 <- plot_grid(p_intercepted, p_indoor, p_exotic, nrow=3)
fig7
save_plot("Figure7.pdf", fig7, base_height=8, base_width=8)
save_plot("Figure7.png", fig7, base_height=8, base_width=8)
```

### PGLS Models
As we did for the legumes, we ran PGLS models to account for  phylogenetic relationships between ant species when measuring the effects of mutualisms on introduction success and introduced range size. The problem with PGLS models is that we lose a lot of trait data for species not in the phylogeny. We either need a better phylogeny or a non-phylogenetic model. 

```{r efns/domatia and ants models, warning=FALSE, message=FALSE}
#all ant mutualisms
#Native range size
ant_tree <- read.tree("Ant_tree.tre") #Reading in ant phylogeny
phy_int <- intersect(ant_tree$tip.label, area$Phy)
phy_diff <- setdiff(ant_tree$tip.label, area$Phy)
pruned_ant_tree <- drop.tip(ant_tree, as.character(phy_diff)) #pruning tree to contain only tips in the dataset

area$n_introduced_ranges <- ifelse(area$total.area.introduced == 0, 0, area$n_introduced_ranges)

area_phy <- area[area$Phy %in% phy_int, c("Phy", "abs_lat_native", "total.area.introduced", "total.area.native", "n_introduced_ranges", "EFN", "Domatia", "Seed_Dispersal")] #dataset for pgls
area_phy <- area_phy[complete.cases(area_phy), ]

#All mutualistic interactions
a11 <- gls(((total.area.native)^(1/4)) ~ EFN + Domatia + Seed_Dispersal + abs_lat_native, 
                correlation = corPagel(1, phy = pruned_ant_tree, form = ~ Phy), 
                method = "ML", data = area_phy) 
summary(a11)

#Numberof non-contiguous ranges
a3 <- gls(log(n_introduced_ranges+1) ~  EFN + Domatia + Seed_Dispersal +  total.area.native + abs_lat_native,  correlation = corPagel(1, phy = pruned_ant_tree, form = ~Phy), method = "ML", data = area_phy)  
summary(a3)

```
#### Pagel correlations between traits
As we did for the legumes, we calculated Pagel's correlation for each pair of mutualistic traits to understand if correlated evolution occurred.

```{r ant Pagel, warning=FALSE, message=FALSE, eval=FALSE}
#Subsetting data to necessary columns
area_pagel <- area[area$Phy %in% area$Phy, c("Phy", "EFN", "Domatia", "Seed_Dispersal")]
area_pagel <- area_pagel[complete.cases(area_pagel), ]

ant_tree <- read.tree("Ant_tree.tre") #Reading in ant phylogeny
phy_int1 <- intersect(ant_tree$tip.label, area_pagel$Phy)
phy_diff1 <- setdiff(ant_tree$tip.label, area_pagel$Phy)
ant_tree_pagel <- drop.tip(ant_tree, as.character(phy_diff1)) #Pruning tree to contain tips in the dataset

area_pagel <- area_pagel[area_pagel$Phy %in% ant_tree_pagel$tip.label,]

ant_traits <- data.frame(area_pagel$Phy, area_pagel$EFN, area_pagel$Domatia, area_pagel$Seed_Dispersal)
colnames(ant_traits) <- c("Phy", "EFN", "Domatia", "Seed_dispersal")
ant_traits <- arrange(ant_traits, ant_tree_pagel$tip.label)


rownames(ant_traits) <- ant_traits[,1]
ant_traits[,1] <- NULL
head(ant_traits)

EFN <-setNames(ant_traits$EFN,rownames(ant_traits))
Domatia <-setNames(ant_traits$Domatia,rownames(ant_traits))
sd <-setNames(ant_traits$Seed_dispersal,rownames(ant_traits))

#EFN and Domatia
fit.aefndom <- fitPagel(ant_tree_pagel, EFN, Domatia) #stat sig
aefndom <- plot(fit.aefndom,lwd.by.rate=TRUE)


fit.adomsd <- fitPagel(ant_tree_pagel, Domatia, sd) #not sig


fit.asdefn <- fitPagel(ant_tree_pagel, sd, EFN) #pval 0.04
asdefn <- plot(fit.asdefn,lwd.by.rate=TRUE)


```

### Latitudinal distribution
We visualized the latitudinal distribution of each mutualism to understand how different mutualisms are distributed globally. For traits where we have data for both legumes and ant (e.g. the EFN-ant interaction), we plotted the distribution of both ants and legumes together.

```{r Latitude}
#devtools::install_github("katiejolly/nationalparkcolors")

library(nationalparkcolors)
pal="Saguaro"
color_EFN = park_palette(pal)[1]
color_dom = park_palette(pal)[2]
color_seed = park_palette(pal)[3]
color_fixer = park_palette(pal)[4]
color_AM = park_palette(pal)[5]
color_EM = park_palette(pal)[6]
x_limit = c(0, 90)
y_limit = c(-0.1, 0.1)
a = 0.5
x_label = 0.05
y_label = 75

#Models
lmer_lat_ant <- lmer(abs_lat_native~EFN+Domatia+Seed_Dispersal+introducedY+(1|tribe), data=area)
summary(lmer_lat_ant)
Anova(lmer_lat_ant, type=3)

lmer_lat_plant <- lmer(abs_lat_native~EFN+Domatia+fixer+introducedY+(1|tribe), data=df)
summary(lmer_lat_plant)
Anova(lmer_lat_plant, type=3)

#Figures
  
#EFN
p_lat_EFN <- ggplot()+
  geom_density(data=subset(area, introducedY == 1 & EFN == 1), aes(y=abs_lat_native,x=after_stat(density)), fill=color_EFN, color=color_EFN, alpha=a)+
  geom_density(data=subset(df, introducedY == 1 & EFN == 1), aes(y=abs_lat_native, x=-after_stat(density)), fill=color_EFN, color=color_EFN, alpha=a)+
  geom_segment(data=subset(area, introducedY == 1 & EFN == 1), aes(y=median(abs_lat_native, na.rm=TRUE), yend=median(abs_lat_native, na.rm=TRUE), x=0, xend=Inf), color=color_EFN, linetype="dashed")+
  geom_segment(data=subset(df, introducedY == 1 & EFN == 1), aes(y=median(abs_lat_native, na.rm=TRUE), yend=median(abs_lat_native, na.rm=TRUE), x=0, xend=-Inf), color=color_EFN, linetype="dashed")+
  annotate("text", x=x_label, y=y_label, label="EFN-visiting ants")+
  annotate("text",x=-x_label, y=y_label, label="EFN-bearing legumes")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_EFN

#Domatia
p_lat_dom <- ggplot()+
  geom_density(data=subset(area, introducedY == 1 & Domatia == 1), aes(y=abs_lat_native,x=after_stat(density)), fill=color_dom, color=color_dom, alpha=a)+
  geom_density(data=subset(df, introducedY == 1 & Domatia == 1), aes(y=abs_lat_native, x=-after_stat(density)), fill=color_dom, color=color_dom, alpha=a)+
  geom_segment(data=subset(area, introducedY == 1 & Domatia == 1), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=Inf), color=color_dom, linetype="dashed")+
  geom_segment(data=subset(df, introducedY == 1 & Domatia == 1), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=-Inf), color=color_dom, linetype="dashed")+
  annotate("text", y=y_label, x=x_label, label="Domatia-nesting ants")+
  annotate("text", y=y_label, x=-x_label, label="Domatia-bearing legumes")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_dom

#Myrmecochory
p_lat_seed <- ggplot()+
  geom_density(data=subset(area, introducedY == 1 & Seed_Dispersal == 1), aes(y=abs_lat_native,x=after_stat(density)), fill=color_seed, color=color_seed, alpha=a)+
  geom_segment(data=subset(area, introducedY == 1 & Seed_Dispersal == 1), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=Inf), color=color_seed, linetype="dashed")+
  annotate("text", y=y_label, x=x_label, label="Seed-dispersing ants")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_seed

#Nodules
p_lat_fixer <- ggplot()+
  geom_density(data=subset(df, introducedY == 1 & fixer == 1), aes(y=abs_lat_native, x=-after_stat(density)), fill=color_fixer, color=color_fixer, alpha=a)+
  geom_segment(data=subset(df, introducedY == 1 & fixer == 1), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=-Inf), color=color_fixer, linetype="dashed")+
  annotate("text", y=y_label, x=-x_label, label="Nodulating legumes")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_fixer

#AM fungi
p_lat_AM <- ggplot()+
  geom_density(data=subset(df, introducedY == 1 & AM == "Y"), aes(y=abs_lat_native, x=-after_stat(density)), fill=color_AM, color=color_AM, alpha=a)+
  geom_segment(data=subset(df, introducedY == 1 & AM == "Y"), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=-Inf), color=color_AM, linetype="dashed")+
  annotate("text", y=y_label, x=-x_label, label="AM legumes")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_AM

#EM fungi
p_lat_EM <- ggplot()+
  geom_density(data=subset(df, introducedY == 1 & EM == "Y"), aes(y=abs_lat_native, x=-after_stat(density)), fill=color_EM, color=color_EM, alpha=a)+
  geom_segment(data=subset(df, introducedY == 1 & EM == "Y"), aes(y=mean(abs_lat_native, na.rm=TRUE), yend=mean(abs_lat_native, na.rm=TRUE), x=0, xend=-Inf), color=color_EM, linetype="dashed")+
  annotate("text", y=y_label, x=-x_label, label="EM legumes")+
  geom_vline(aes(xintercept=0))+
  scale_x_continuous(limits=y_limit)+
  scale_y_continuous(limits=x_limit)+
ylab(expression(paste("Latitude ", "|", "\u00B0", "|")))+
  xlab("Density")+
  theme_cowplot()
p_lat_EM

p_lat <- plot_grid(p_lat_EFN, p_lat_dom, p_lat_seed, p_lat_fixer, p_lat_AM, p_lat_EM, nrow=6, labels="AUTO")
p_lat
save_plot("Latitude Figure.pdf", p_lat, base_height = 12, base_width=6)
```