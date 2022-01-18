#install.packages("taxize")
library(taxize)

#Finding family name of soud dataset

soud <- read.csv("Soudzilovskaia_analysis.csv")
mahe <- read.csv("Mycorrhizae_species_list.csv")

#dat.fam <- tax_name(sci= soud$Phy, get = 'family', db = 'ncbi')


#final <- data.frame(comb, dat.tribe$tribe)
#write.csv(final, "combined396_withtribes.csv")

range <- read.csv("legume_invasion_data_simonsenetal.csv")

soud_subs <- soud[soud$Phy %in% range$Species, ] #subsetting soud by species that are present in range
write.csv(soud_subs, "Soudzilovskaia_simonsen_intersect.csv")
mahe_subs <- mahe[mahe$ï..phy %in% range$Species, ] #subsetting mahe by species that are present in range
write.csv(mahe_subs, "Maherali_simonsen_intersect.csv")
