library(viridis)
library(data.table)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(zoo)
library(mgcv)
library(scam)
library(scales)


dir.create('FIGURES')
dir.create('NUMBERS')

theme_common <- theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = 14),
    axis.title       = element_text(size = 16),
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 14)
  )



file_path <- "COL_2026-02-13_XR/NameUsage.tsv"


nameusage <- fread(file_path, sep = "\t", quote = "", data.table = FALSE)
accepted_spp<-nameusage[which(nameusage$`col:status`=='accepted'),]
animals<-accepted_spp[which(accepted_spp$`col:kingdom`=='Animalia'),]
animals<-animals[animals$`col:rank`=='species',]
n_all<-nrow(animals)
animals <- animals[which(animals$`col:extinct` == FALSE | is.na(animals$`col:extinct`)), ]
n_extant<-nrow(animals)
n_extinct<-n_all-n_extant

#identify year of first description (keeping into account potential taxonomic revisions)
animals$`col:basionymAuthorshipYear`[animals$`col:basionymAuthorshipYear`=='']<-NA
animals$`col:publishedInYear`[animals$`col:publishedInYear`=='']<-NA
animals$`col:combinationAuthorshipYear`[animals$`col:combinationAuthorshipYear`=='']<-NA

animals$`col:basionymAuthorshipYear`<-as.numeric(animals$`col:basionymAuthorshipYear`)
animals$`col:publishedInYear`<-as.numeric(animals$`col:publishedInYear`)
animals$`col:combinationAuthorshipYear`<-as.numeric(animals$`col:combinationAuthorshipYear`)

first_descr_animals<-coalesce(animals$`col:basionymAuthorshipYear`,animals$`col:publishedInYear`,animals$`col:combinationAuthorshipYear`)
first_descr_animals <- substr(first_descr_animals, 1, 4)

#make more readable table with only relevant columns
animals<-animals[which(!is.na(first_descr_animals)),]
animals<-cbind(animals$`col:scientificName`,first_descr_animals[which(!is.na(first_descr_animals))],
               animals$`col:order`,
               animals$`col:subclass`,
               animals$`col:class`,
               animals$`col:subphylum`,
               animals$`col:phylum`
               )



animals<-data.frame(animals)
colnames(animals)<-c('species','year','order','subclass','class','subphylum','phylum')

animals$year<-as.numeric(animals$year)

#keep only descriptions from 1500 on
animals<-animals[which(animals$year>1500),]

#remove obvious errors (species with future year of description)
animals<-animals[which(animals$year<2027),]
nrow(animals)



####ESTIMATE TRENDS
sp_per_y<-aggregate(animals$species~animals$year,FUN='length')
sp_per_y<-sp_per_y[which(sp_per_y$`animals$year`>=1970),]
sp_per_y<-sp_per_y[which(sp_per_y$`animals$year`<=2020),]
df<-data.frame('year'=as.numeric(as.character(sp_per_y$`animals$year`)),
               'descriptions'=sp_per_y$`animals$species`)


df$year<-df$year-min(df$year)
lin_mod_rec_desc<-lm(df$descriptions~df$year)
summary(lin_mod_rec_desc)

pdf('./FIGURES/description_trend.pdf',width=4.5,height=4.5)
plot(df$year+1970,df$descriptions/1000,
     xlab='year',ylab='descriptions (k species)',las=1,cex.axis=1.2,cex.lab=1.2,
     pch=16)

y_adj<-df$year+1970
abline(lm(df$descriptions/1000~y_adj))
dev.off()



a<-lin_mod_rec_desc$coefficients[1]
b<-lin_mod_rec_desc$coefficients[2]

a+(2026-1970)*b
sink('./NUMBERS/description_rate_equation.txt')
print(summary(lin_mod_rec_desc))
sink()



##################################################
##################################################
###Extinction rates

iucn_dir <- './iucn_species_data'
folders <- list.dirs(iucn_dir, recursive = FALSE)

# Read and combine all assessments and taxonomy files
assessments <- rbindlist(lapply(folders, function(f) fread(file.path(f, "assessments.csv"))))
taxonomy    <- rbindlist(lapply(folders, function(f) fread(file.path(f, "taxonomy.csv"))))

# Remove duplicates in case species appear in multiple folders
assessments <- unique(assessments, by = "internalTaxonId")
taxonomy    <- unique(taxonomy,    by = "internalTaxonId")

# Join, keep only animals, exclude DD
iucn_sp <- merge(assessments, taxonomy, by = "internalTaxonId")
iucn_sp <- iucn_sp[iucn_sp$kingdomName == "ANIMALIA", ]
iucn_eval <- iucn_sp[iucn_sp$redlistCategory != "Data Deficient", ]


#use this code to identify which "yearLastSeen" values need to be fixed/converted to year
# year_last_seen<-unique(iucn_eval$yearLastSeen)
# write.table(unique(year_last_seen[is.na(as.numeric(year_last_seen))]),'y_last_seen_to_fix.txt',
#             row.names=F,col.names=F,quote=F)


####the following dictionary was obtained by looking at uniques entries in extinct_species last seen dates
last_seen_to_year <- c(
  "early 1900s"                                              = 1910,
  "1900s"                                                    = 1909,
  "post 1950s"                                               = 1960,
  "1970s"                                                    = 1979,
  "August 1975"                                              = 1975,
  "known only from recent fossil deposits"                   = -9999,
  "presumably following European arrival"                    = 1500,  # too vague
  "1935 (confirmed); 1975 (unconfirmed)"                     = 1975,
  "1870s"                                                    = 1879,
  "1945-1947"                                                = 1947,
  "last captive 1938; last wild animal known in 1932"        = 1938,
  "1860-1862"                                                = 1862,
  "1960s"                                                    = 1969,
  "unknown"                                                  = 1500,
  "early 18th century"                                       = 1720,
  "March 17, 1984"                                           = 1984,
  "end of December 1863"                                     = 1863,
  "1940?"                                                    = 1940,
  "19th century?"                                            = 1900,
  "20th century?"                                            = 2000,
  "1950s"                                                    = 1959,
  "ca 1850"                                                  = 1850,
  "post 1927"                                                = 1927,
  "pre 1911"                                                 = 1911,
  "1799-1800"                                                = 1800,
  "mid-19th Century"                                         = 1850,
  "extinction estimated in year 1525"                        = 1525,
  "May 1967"                                                 = 1967,
  "19th Century"                                             = 1900,
  "N/A"                                                      = 1500,
  "1870-1878"                                                = 1878,
  "Unknown: possibly 1950s or 1970."                         = 1970,
  "Never recorded live"                                      = -9999,
  "Fossil deposit"                                           = -9999,
  "Fossil deposits"                                          = -9999,
  "1930s"                                                    = 1939,
  "1901-1902"                                                = 1902,
  "1843-1844"                                                = 1844,
  "1840s"                                                    = 1849,
  "early 1950s (confirmed)"                                  = 1955,
  "pre 1900"                                                 = 1900,
  "1980s"                                                    = 1989,
  "1990s"                                                    = 1999,
  "1943 (confirmed); 1960s (unconfirmed)"                    = 1969,
  "early 20th century"                                       = 1920,
  "1700-1800"                                                = 1800,
  "1800-1878"                                                = 1878,
  "1874 (with certainty); observations from the 1890s"       = 1899,
  "before 1874"                                              = 1874,
  "1968 (confirmed); possibly seen in late 1970s"            = 1979,
  "1897-1898"                                                = 1898,
  "March 1985"                                               = 1985,
  "Known only from recent fossil deposits"                   = -9999,
  "1933 (confirmed)"                                         = 1933,
  "1890s"                                                    = 1899,
  "1930?"                                                    = 1930,
  "Possibly 1803"                                            = 1803,
  "uncertain"                                                = 1500,
  "known only from fossils"                                  = -9999,
  "before 1856"                                              = 1856,
  "before 1947"                                              = 1947,
  "July 1862"                                                = 1862,
  "1989/90"                                                  = 1990,
  "before 1904"                                              = 1904,
  "Late 1800s"                                               = 1899,
  "1940s-1980s"                                              = 1989,
  "2016 (questionable)"                                      = 2016,
  "1920s"                                                    = 1929,
  "Early 1900s"                                              = 1910,
  "before 1859"                                              = 1859,
  "Fossil"                                                   = -9999,
  "pre 1658"                                                 = 1658,
  "around 1700"                                              = 1700,
  "Subfossil record"                                         = -9999,
  "27 August 2009"                                           = 2009,
  "There is no record of it ever having been seen."          = -9999,
  "before 1960"                                              = 1960,
  "1450 BP (estimated year of extinction)"                   = -9999,  # ~500 AD, treated as fossil
  "1933-1935"                                                = 1935,
  "July 1878"                                                = 1878,
  "1904-1908"                                                = 1908,
  "19th century"                                             = 1900,
  "post-1500"                                                = 1500,
  "post-1550"                                                = 1550,
  "late 20th century"                                        = 1999,
  "1825-1854"                                                = 1854,
  "1827-1828"                                                = 1828,
  "1726-1761"                                                = 1761,
  "1895-1900"                                                = 1900,
  "1850-1900"                                                = 1900,
  "end of 19th century, probably 1871"                       = 1899,
  "never recorded alive"                                     = -9999,
  "1863-1873"                                                = 1873,
  "1600s"                                                    = 1609,
  "1250-1860"                                                = 1860,
  "1600-1700"                                                = 1700,
  "17th Century"                                             = 1700,
  "Quaternary"                                               = -9999,  # geological period
  "1850s"                                                    = 1859,
  "1464\u20131636"                                           = 1636,  # en-dash variant
  "1839-1841"                                                = 1841,
  "1908-1909"                                                = 1909,
  "11/04/1890"                                               = 1890,
  "26/08/1967"                                               = 1967,
  "late nineteenth century"                                  = 1899,
  "Early 20th century"                                       = 1920,
  "Mid 1900s"                                                = 1950,
  "2003 (dead)"                                             = 2003,
  "1780-1800"                                               = 1800,
  "Unknown"                                                 = 1500,
  "pre 1889"                                                = 1889,
  "last absolute confirmed date is 1969/70"                 = 1970,
  "before 1998"                                             = 1998,
  "1860s"                                                   = 1869,
  "2000-2006"                                               = 2006,
  "before 1992"                                             = 1992,
  "(1934)"                                                  = 1934,
  "around 1970"                                             = 1970,
  "1860-1867"                                               = 1867,
  "1850-1859 (dead)"                                        = 1859,
  "(1880s)"                                                 = 1889,
  "late 1980s-early 1990s"                                  = 1993,
  "1968 (confirmed); 2016 (unconfirmed)"                    = 2016,
  "31/08/1985"                                              = 1985,
  "1886-1888 (confirmed); 1960s (possible)"                 = 1969,
  "presumably 1886-1888"                                    = 1888,
  "none certain"                                            = 1500,
  "19th or early 20th Century"                              = 1920,
  "before 1893"                                             = 1893,
  "before 1960s"                                            = 1960,
  "June 10, 1968"                                           = 1968,
  "1986-1988"                                               = 1988,
  "2001-2005"                                               = 2005,
  "mid-1980s"                                               = 1985,
  "mid 1980s"                                               = 1985,
  "early 1990s"                                             = 1993,
  "Before 1956"                                             = 1956,
  "Before 2009"                                             = 2009,
  "9 August 1976"                                           = 1976,
  "24 December 1940"                                        = 1940,
  "10 August 1964"                                          = 1964,
  "March 2001"                                              = 2001,
  "before 2000"                                             = 2000,
  "before 1973"                                             = 1973,
  "early 1980s"                                             = 1983,
  "late 1970s"                                              = 1979,
  "dead shells between 2000 and 2011"                       = 2011,
  "late 1990s"                                              = 1999,
  "never seen in the wild"                                  = -9999,
  "before 1983"                                             = 1983,
  "before 2003"                                             = 2003,
  "5 July 1976"                                             = 1976,
  "2000s"                                                   = 2009,
  "before 1994"                                             = 1994,
  "1968-1985"                                               = 1985,
  "1830s"                                                   = 1839,
  "before 1879"                                             = 1879,
  "1857-1858"                                               = 1858,
  "28 February 1941"                                        = 1941,
  "22 September 2005"                                       = 2005,
  "March 1967"                                              = 1967,
  "February 1876"                                           = 1876,
  "05 March 1967"                                           = 1967,
  "12-19 March 1967"                                        = 1967,
  "26 November 1965"                                        = 1965,
  "December 1965"                                           = 1965,
  "01 December 1965"                                        = 1965,
  "1/4/1967"                                                = 1967,
  "late 19th century"                                       = 1899,
  "2010 (but not living)"                                   = 2010,
  "2012 dead shells"                                        = 2012,
  "ca 1900"                                                 = 1900,
  "1920-1929"                                               = 1929,
  "July 2010"                                               = 2010,
  "2000-2013"                                               = 2013,
  "before 1967"                                             = 1967,
  "06/1928"                                                 = 1928,
  "March 27th 2006"                                         = 2006,
  "22 June 2010"                                            = 2010,
  "June 27 - 2010"                                          = 2010,
  "5 August 1988"                                           = 1988,
  "2010s"                                                   = 2019
)



convert_year <- function(x, dict) {
  result <- dict[x]                                        # lookup dizionario
  not_found <- is.na(result)
  result[not_found] <- suppressWarnings(as.numeric(x[not_found]))  # anni come stringa ("1985" ecc.)
  result[is.na(result) | x == ""] <- 1500                 # NA residui ed empty → 1500
  result
}


###EVALUATE SPECIES PER YEAR
iucn_eval$year_num <- convert_year(iucn_eval$yearLastSeen, last_seen_to_year)

###then exclude species from categories EX, EW, and CR+PE with last_seen<1900 (-9999 = fossil record)
extinct_species<-iucn_eval[iucn_eval$redlistCategory %in% c("Extinct", "Extinct in the Wild"),]$scientificName.x
critically_endangered_species<-iucn_eval[iucn_eval$redlistCategory=="Critically Endangered",]$scientificName.x
possibly_extinct_species<-iucn_eval[iucn_eval$possiblyExtinct == TRUE,]$scientificName.x
last_seen_pre1900<-iucn_eval[iucn_eval$year_num<1900,]$scientificName.x

to_exclude<-intersect(c(extinct_species,intersect(critically_endangered_species,
                                                  possibly_extinct_species)),last_seen_pre1900)


iucn_eval<-iucn_eval[!(iucn_eval$scientificName.x %in% to_exclude),]

n_eval <- nrow(iucn_eval)

extinct_species<-iucn_eval[iucn_eval$redlistCategory %in% c("Extinct", "Extinct in the Wild"),]$scientificName.x
critically_endangered_species<-iucn_eval[iucn_eval$redlistCategory=="Critically Endangered",]$scientificName.x
possibly_extinct_species<-iucn_eval[iucn_eval$possiblyExtinct == TRUE,]$scientificName.x



n_ext_l<-length(extinct_species)
n_ext_m <-n_ext_l+length(intersect(critically_endangered_species,possibly_extinct_species)) 
n_ext_h <-n_ext_m+length(critically_endangered_species)


ext_rate_l <- (n_ext_l / n_eval) / (2025-1900)
ext_rate_m <- (n_ext_m / n_eval) / (2025-1900)
ext_rate_h <- (n_ext_h / n_eval) / (2025-1900)


background_ER <- 2 / 10^6
sink('./NUMBERS/iucn_ext_rate_vs_BG.txt')
cat("IUCN Rate to Background Ratio:\n")
cat("Low Scenario:   ", ext_rate_l / background_ER, "\n")
cat("Middle Scenario:", ext_rate_m / background_ER, "\n")
cat("High Scenario:  ", ext_rate_h / background_ER, "\n")
sink()


####extinction rates by taxon (class) for Fig. 1
# ── Extinction rates by class ─────────────────────────────────────────────────

classes <- unique(iucn_eval$className)

ext_rate_by_class <- rbindlist(lapply(classes, function(cl) {
  
  sub <- iucn_eval[iucn_eval$className == cl, ]
  
  n_eval_cl <- nrow(sub)
  
  ext_sp_cl     <- sub[sub$redlistCategory %in% c("Extinct", "Extinct in the Wild"), ]$scientificName.x
  cr_sp_cl      <- sub[sub$redlistCategory == "Critically Endangered", ]$scientificName.x
  pe_sp_cl      <- sub[sub$possiblyExtinct == TRUE, ]$scientificName.x
  
  n_ext_l_cl <- length(ext_sp_cl)
  n_ext_m_cl <- n_ext_l_cl + length(intersect(cr_sp_cl, pe_sp_cl))
  n_ext_h_cl <- n_ext_m_cl + length(cr_sp_cl)
  
  data.table(
    class        = cl,
    n_eval       = n_eval_cl,
    n_extinct_l  = n_ext_l_cl,
    n_extinct_m  = n_ext_m_cl,
    n_extinct_h  = n_ext_h_cl,
    ext_rate_l   = (n_ext_l_cl / n_eval_cl) / (2025 - 1900),
    ext_rate_m   = (n_ext_m_cl / n_eval_cl) / (2025 - 1900),
    ext_rate_h   = (n_ext_h_cl / n_eval_cl) / (2025 - 1900)
  )
}))

# Sort by number of evaluated species
ext_rate_by_class <- ext_rate_by_class[order(-n_eval)]

# Save
dir.create('./NUMBERS', showWarnings = FALSE)
fwrite(ext_rate_by_class, './NUMBERS/ext_rate_by_class.csv')



##############################################################
###Time to Mass Extinction under different extinction rates

lambda <- 1e-7 ###Speciation rate
ext_rate_range<-range(ext_rate_l,ext_rate_h)
mu_intervals <- seq(ext_rate_range[1], ext_rate_range[2], length.out = 10)

plot_points <- data.frame(
  extinction_rate = mu_intervals,
  years_to_75 = log(0.25) / (lambda - mu_intervals)
) 


pA <- ggplot(plot_points, aes(x = years_to_75, y = extinction_rate,
                              color = factor(extinction_rate))) +
  geom_point(size = 4, alpha = 0.9) +
  geom_line(aes(group = 1), linetype = "solid", alpha = 0.5) +
  scale_color_viridis_d(
    name   = "extinction rate",
    option = "viridis",
    labels = function(x) formatC(as.numeric(x), format = "e", digits = 2),
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_x_log10(breaks = c(2500, 5000, 10000, 20000, 40000)) +
  labs(
    x = "years to 75% loss",
    y = "extinction rate"
  ) +
  
  theme_bw() +
  theme_common +
  theme(
    legend.position   = c(0.8, 0.55),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )


t_max <- max(log(0.01) / (lambda - mu_intervals))
t_seq <- seq(1, t_max, length.out = 1000)

curves_df <- lapply(mu_intervals, function(mu) {
  r <- lambda - mu
  data.frame(
    year            = t_seq + 2026,
    diversity_pct   = pmax(100 * exp(r * t_seq), 0),
    extinction_rate = mu,
    mu_factor       = factor(mu)
  )
}) %>% bind_rows()

arrival_pts <- plot_points %>%
  mutate(
    year          = years_to_75 + 2026,
    diversity_pct = 25,
    mu_factor     = factor(extinction_rate) 
  )

pB <- ggplot(curves_df,
             aes(x = year, y = diversity_pct,
                 group = mu_factor,
                 color = extinction_rate)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 25, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +

  geom_point(data = arrival_pts,
             aes(x = year, y = diversity_pct,
                 group = mu_factor),
             size = 3.5, alpha = 0.9) +
  scale_color_viridis_c(name = "extinction rate",
                        option = "viridis") +
  scale_x_continuous(
    trans  = "log10",
    breaks = c(2100, 5000, 10000, 20000, 40000,80000),
    labels = c("2100", "5000", "10000", "20000", "40000", "80000")
  ) +
  scale_y_continuous(limits = c(0, 105),
                     breaks = seq(0, 100, 25),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x = "year",
    y = "remaining diversity (%)"
  ) +
  theme_bw() +
  theme_common +
  theme(
    legend.position   = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    plot.margin       = margin(5, 40, 5, 5)
  )+
  coord_cartesian(clip = "off")


pdf("./FIGURES/time_to_loss_paired.pdf",width=12,height=5)

print(
  pA + pB +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag          = element_text(size = 18, face = "bold"),
      plot.tag.position = "top"   # centrato sopra il pannello
    )
)
dev.off()




sink('./NUMBERS/time_to_ME_summary.txt')
summary(plot_points$years_to_75)+2026
sink()


##########################################################
#Global diversity ranges
costello<-c(1800000,2000000,1900000)
larsen<-c(15300000,163200000,89250000)

div_range<-range(c(costello,larsen))

#################
get_descriptions_y <- function(year, a, b, C_prev, S, alpha=0.5, base_year=1970) {
  r_linear <- a + b * (year - base_year)
  depletion <- max(0, 1 - C_prev / S)
  depletion_mod <- depletion ^ alpha
  pred <- max(0, r_linear * depletion_mod)
  return(pred)
}

y_0 <- 2025
desc_n_y0 <- sum(animals$year <= y_0)
y_end <- 5000
tot <- 7770000  # Mora 2011
alpha_fixed <- 0.5

# --- Plot 1: fixed S, varying alpha ---
alpha_vals <- seq(0.0, 1, 0.01)
desc_traj <- c()
for (alpha in alpha_vals) {
  tot_des <- desc_n_y0
  de <- c()
  de_y <- c()
  for (i in (y_0:y_end)) {
    des_y <- get_descriptions_y(i - y_0, a, b, tot_des, tot, alpha)
    if (des_y >= (tot - tot_des)) { des_y <- (tot - tot_des) }
    de_y <- c(de_y, des_y)
    tot_des <- tot_des + des_y
    de <- c(de, tot_des)
  }
  desc_traj <- rbind(desc_traj, cbind(alpha, y_0:y_end, de_y, de/tot))
  print(alpha)
}


div_range<-range(c(costello,larsen))
div_range<-seq(min(div_range,desc_n_y0),max(div_range),length.out=100)
# --- Plot 2: fixed alpha=0.5, varying S ---
div_range<-c(div_range,7.77e6)
div_range<-sort(div_range)
s_traj <- c()
for (S_val in div_range) {
  tot_des <- desc_n_y0
  de <- c()
  de_y <- c()
  for (i in (y_0:y_end)) {
    des_y <- get_descriptions_y(i - y_0, a, b, tot_des, S_val, alpha_fixed)
    if (des_y >= (S_val - tot_des)) { des_y <- max(0, S_val - tot_des) }
    de_y <- c(de_y, des_y)
    tot_des <- tot_des + des_y
    de <- c(de, tot_des)
  }
  s_traj <- rbind(s_traj, cbind(S_val, y_0:y_end, de_y, de/S_val))
  print(S_val)
}




alpha_df <- do.call(rbind, lapply(seq_along(alpha_vals), function(m) {
  alpha <- alpha_vals[m]
  d <- desc_traj[desc_traj[, 1] == alpha, ]
  data.frame(
    x     = 100 * d[, 4],
    y     = 100 * d[, 3] / tot,
    alpha = alpha,
    highlight = FALSE
  )
}))

alpha_fixed_df <- {
  d <- desc_traj[desc_traj[, 1] == alpha_fixed, ]
  data.frame(x = 100 * d[, 4], y = 100 * d[, 3] / tot,
             alpha = alpha_fixed, highlight = TRUE)
}

pa <- ggplot(alpha_df, aes(x = x, y = y, group = alpha, color = alpha)) +
  geom_line() +
  geom_line(data = alpha_fixed_df,
            aes(x = x, y = y),
            color = "red", linewidth = 1.2, linetype = "dashed",
            inherit.aes = FALSE) +
  scale_color_viridis_c(
    name   = expression(alpha),
    breaks = seq(0, 1, 0.2)
  ) +
  labs(x = "% described global diversity",
       y = "% global diversity described per year") +
  theme_common+
  theme(
    legend.position   = c(0.01, 0.99),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

s_df <- do.call(rbind, lapply(seq_along(div_range), function(m) {
  S_val <- div_range[m]
  d <- s_traj[s_traj[, 1] == S_val, ]
  data.frame(
    x     = 100 * d[, 4],
    y     = 100 * d[, 3] / S_val,
    S_val = S_val
  )
}))

mora_df <- {
  d <- s_traj[s_traj[, 1] == 7.77e6, ]
  data.frame(x = 100 * d[, 4], y = 100 * d[, 3] / 7.77e6)
}

pb <- ggplot(s_df, aes(x = x, y = y, group = S_val, color = S_val)) +
  geom_line(linewidth = 1) +
  geom_line(data = mora_df,
            aes(x = x, y = y),
            color = "red", linewidth = 1.2, linetype = "dashed",
            inherit.aes = FALSE) +
  scale_color_viridis_c(
    name   = "global diversity\n(million species)",
    breaks = div_range[round(seq(1, length(div_range), length.out = 5))],
    labels = function(x) format(x / 1e6, digits = 2)
  ) +
  labs(x = "% described global diversity",
       y = "% global diversity described per year") +
  theme_common+
  theme(
    legend.position   = c(0.01, 0.99),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

# --- Combina e salva ---
pdf("./FIGURES/alpha_vals.pdf", width = 12, height = 5)
print(
  pa + pb +
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag          = element_text(size = 18, face = "bold"),
          plot.tag.position = "top")
)
dev.off()



#########################################################
####Complete model with extinctions and descriptions
spec_r <- 1e-07
y_0 <- 2026
y_end <- 100000
desc_n_y0 <- sum(animals$year < y_0)
div_range <- range(c(costello, larsen))

res_trajectories <- list()
res_thresholds   <- list()

ext_grid <- seq(ext_rate_range[1], ext_rate_range[2], length.out = 10)
tot_grid <- seq(div_range[1],      div_range[2],      length.out = 10)

for (ext_r in ext_grid) {
  for (tot in tot_grid) {
    
    n_years <- y_end - y_0 + 1
    div <- numeric(n_years); des <- numeric(n_years)
    ue  <- numeric(n_years); doc_ext <- numeric(n_years)
    des_alive <- numeric(n_years)
    
    y_div25 <- NA; y_doc75 <- NA; y_ue25 <- NA
    y_des75 <- NA
    div[1]       <- tot
    des[1]       <- desc_n_y0
    des_alive[1] <- desc_n_y0
    ue[1]        <- 0
    doc_ext[1]   <- 0
    
    for (t in 2:n_years) {
      prev <- t - 1
      
      div[t]     <- div[prev] * (1 + spec_r - ext_r)
      des_y      <- get_descriptions_y(year = y_0 + t - 1, a = a, b = b,
                                       C_prev = des[prev], S = tot, alpha = 0.5)
      des_y      <- max(0, min(des_y, div[prev] - des[prev]))
      doc_ext[t] <- doc_ext[prev] + des_alive[prev] * ext_r
      des_alive[t] <- des_alive[prev] * (1 - ext_r) + des_y
      des[t]     <- des[prev] + des_y
      ue[t]      <- ue[prev] + (div[prev] - des_alive[prev]) * ext_r
      
      if (is.na(y_div25)      && div[t] < 0.25 * tot)       y_div25      <- y_0 + t - 1
      if (is.na(y_doc75)      && doc_ext[t] > 0.75 * tot)   y_doc75      <- y_0 + t - 1
      if (is.na(y_ue25)       && ue[t] > 0.25 * tot)        y_ue25       <- y_0 + t - 1
      if (is.na(y_des75)      && des[t] > 0.75 * tot)       y_des75      <- y_0 + t - 1
      all_found <- !any(is.na(c(y_div25, y_doc75, y_ue25, y_des75)))
      if (t > (10000 - y_0) && all_found) break
    }
    
    plot_limit <- min(t, 10000 - y_0 + 1)
    res_trajectories[[length(res_trajectories) + 1]] <- data.frame(
      ext_r = ext_r, tot = tot, year = y_0:(y_0 + plot_limit - 1),
      div   = div[1:plot_limit], des = des[1:plot_limit],
      ue    = ue[1:plot_limit],  doc_ext = doc_ext[1:plot_limit]
    )
    
    res_thresholds[[length(res_thresholds) + 1]] <- data.frame(
      ext_r = ext_r, tot = tot,
      y_div25 = y_div25, y_doc75 = y_doc75,
      y_ue25  = y_ue25,  y_des75 = y_des75,
      max_ue_frac      = max(ue) / tot
    )
  }
}

res        <- bind_rows(res_trajectories)
results_df <- bind_rows(res_thresholds)



summarise_threshold <- function(df, col, n_total) {
  x <- df[[col]]
  x <- x[!is.na(x)]
  data.frame(
    threshold   = col,
    n_reached   = length(x),
    pct_reached = round(length(x) / n_total * 100, 1),
    mean        = round(mean(x)),
    median      = round(median(x)),
    ci_low      = round(quantile(x, 0.025)),
    ci_high     = round(quantile(x, 0.975)),
    min         = round(min(x)),
    max         = round(max(x))
  )
}

n_total <- nrow(results_df)

year_summary <- bind_rows(
  summarise_threshold(results_df, "y_ue25",       n_total),
  summarise_threshold(results_df, "y_des75",      n_total),
  summarise_threshold(results_df, "y_div25",      n_total),
)


write.csv(year_summary,       "./NUMBERS/summary_thresholds.csv",    row.names = FALSE)




sink('./NUMBERS/fraction_simulations_where_undocumented_extinctions_exceeded_25.txt')
sum(results_df$max_ue_frac>0.25)/100
sink()


p1 <- ggplot(results_df, aes(x = ext_r, y = tot, fill = max_ue_frac)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "plasma",
    guide  = guide_colorbar(
      title      = "fraction of undocumented extinctions",
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(6, "cm")
    )
  ) +
  labs(x = "extinction rate", y = "global diversity") +
  theme_common +
  theme(legend.position = "top")

p2 <- ggplot(results_df, aes(x = ext_r, y = tot, fill = y_des75)) +
  geom_tile() +
  scale_fill_viridis_c(
    option   = "plasma",
    na.value = "grey80",
    guide    = guide_colorbar(
      title          = "year described diversity > 75%",
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(6, "cm")
    )
  ) +
  labs(x = "extinction rate", y = "global diversity") +
  theme_common +
  theme(legend.position = "top")

combined_plot <- p1 + p2 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.tag.position = "top")

pdf("./FIGURES/undocumented_extinctions_colorplot.pdf",
    width = 10, height = 5.5)
print(combined_plot)
dev.off()

########FUTURE TRAJECTORIES
res_summary <- res %>%
  mutate(
    frac_div = 1 - div / tot,
    frac_ue = ue / tot,
    frac_des = des / tot
  ) %>%
  group_by(year) %>%
  summarise(across(
    c(frac_div, frac_des, frac_ue),
    list(
      mean = ~mean(., na.rm = TRUE),
      low  = ~mean(.) - 1.96 * (sd(.) / sqrt(n())),
      high = ~mean(.) + 1.96 * (sd(.) / sqrt(n())),
      min  = ~min(., na.rm = TRUE),
      max  = ~max(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))

colnames(res_summary) <- gsub("frac_", "", colnames(res_summary))

res_long <- res_summary %>%
  pivot_longer(
    cols = -year,
    names_to = c("var", ".value"),
    names_pattern = "(.*)_(mean|low|high|min|max)"
  )

# Viridis colors for 3 variables
var_labels <- c(
  "div" = "extinctions",
  "des" = "descriptions",
  "ue"  = "undocumented ext"
)

viridis_cols <- viridis(8)
col_div <- viridis_cols[1]
col_des <- viridis_cols[4]
col_ue  <- viridis_cols[7]
transparency <- 0.4

pa <- ggplot() +
  geom_ribbon(data = res_long[res_long$var == "div", ],
              aes(x = year, ymin = low, ymax = high),
              fill = adjustcolor(col_div, alpha.f = transparency)) +
  geom_ribbon(data = res_long[res_long$var == "des", ],
              aes(x = year, ymin = low, ymax = high),
              fill = adjustcolor(col_des, alpha.f = transparency)) +
  geom_ribbon(data = res_long[res_long$var == "ue", ],
              aes(x = year, ymin = low, ymax = high),
              fill = adjustcolor(col_ue, alpha.f = transparency)) +
  geom_line(data = res_long[res_long$var == "div", ],
            aes(x = year, y = mean), color = col_div, linewidth = 1.5) +
  geom_line(data = res_long[res_long$var == "des", ],
            aes(x = year, y = mean), color = col_des, linewidth = 1.5) +
  geom_line(data = res_long[res_long$var == "ue", ],
            aes(x = year, y = mean), color = col_ue,  linewidth = 1.5) +
  scale_color_manual(
    name   = NULL,
    values = setNames(c(col_div, col_des, col_ue),
                      c("extinctions", "descriptions", "undocumented extinctions"))
  ) +
  geom_line(data = res_long,
            aes(x = year, y = mean, color = var), linewidth = 0) +
  scale_color_manual(
    name   = NULL,
    values = c(div = col_div, des = col_des, ue = col_ue),
    guide  = guide_legend(override.aes = list(linewidth = 5, size = 4)),
    labels = c("descriptions","extinctions","undocumented extinctions")
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "year", y = "fraction of initial diversity") +
  theme_common +
  theme(legend.position = c(0.7, 0.4),
        legend.background = element_rect(fill = "white", color = "grey80",
                                         linewidth = 0.1))

pb <- ggplot(res_long[res_long$var == "div", ],
             aes(x = year)) +
  geom_ribbon(aes(ymin = min, ymax = max),
              fill = adjustcolor(col_div, alpha.f = transparency)) +
  geom_line(aes(y = mean), color = col_div, linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "year", y = "extinctions") +
  theme_common

pc <- ggplot(res_long[res_long$var == "des", ],
             aes(x = year)) +
  geom_ribbon(aes(ymin = min, ymax = max),
              fill = adjustcolor(col_des, alpha.f = transparency)) +
  geom_line(aes(y = mean), color = col_des, linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "year", y = "descriptions") +
  theme_common

pd <- ggplot(res_long[res_long$var == "ue", ],
             aes(x = year)) +
  geom_ribbon(aes(ymin = min, ymax = max),
              fill = adjustcolor(col_ue, alpha.f = transparency)) +
  geom_line(aes(y = mean), color = col_ue, linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Year", y = "undocumented extinctions") +
  theme_common


pdf("./FIGURES/combined_trajectory_plot.pdf", width = 12, height = 10)
print(
  (pa + pb) / (pc + pd) +
    plot_annotation(tag_levels = "a",
                    tag_prefix  = "",
                    tag_suffix  = "") &
    theme(plot.tag = element_text(size = 18, face = "bold"),plot.tag.position = "top")
)
dev.off()



###ESCALATION
lambda        <- 1e-7
ME_threshold  <- 75           

df <- read.csv('extinction_timeline.csv')
df <- df[order(-df$time), ]

loess_fit    <- loess(all_gen_well_res ~ seq_along(time), data = df, span = 0.75)
df$trend     <- predict(loess_fit)
df$detrended <- df$all_gen_well_res - df$trend + mean(df$all_gen_well_res)

window_size <- 25
df$baseline <- rollapply(df$detrended, width = window_size,
                         FUN = function(x) max(x, na.rm = TRUE),
                         fill = "extend", align = "right")

df$perc_loss <- (df$detrended - df$baseline) / df$baseline * 100
df$perc_loss[df$perc_loss > 0] <- 0

df$is_crisis  <- df$perc_loss < -0.5
df$episode_id <- cumsum(df$is_crisis == FALSE & lag(df$is_crisis, default = FALSE))
df$episode_id[!df$is_crisis] <- NA

predictive_df_raw <- df %>%
  filter(!is.na(episode_id)) %>%
  group_by(episode_id) %>%
  mutate(nadir_in_episode = min(perc_loss)) %>%
  filter(row_number() <= which.min(perc_loss)) %>%
  ungroup() %>%
  mutate(
    obs_loss       = abs(perc_loss),
    final_severity = abs(nadir_in_episode)
  )


episode_starts <- predictive_df_raw %>%
  group_by(episode_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(obs_loss = 0)

predictive_df <- bind_rows(predictive_df_raw, episode_starts) %>%
  arrange(episode_id, obs_loss)

#GAM 
results <- data.frame(observed = seq(0, 90, by = 1)) %>%
  rowwise() %>%
  mutate(
    n_at_risk    = sum(predictive_df$obs_loss >= observed),
    n_reached_ME = sum(predictive_df$obs_loss >= observed &
                         predictive_df$final_severity >= ME_threshold),
    prob_ME      = if_else(n_at_risk > 0, n_reached_ME / n_at_risk, NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    prob_ME = if_else(observed >= ME_threshold & !is.na(prob_ME), 1, prob_ME)
  )

train_dat <- results %>%
  rename(loss = observed, prob = prob_ME) %>%
  filter(!is.na(prob), n_at_risk > 0)

gam_mod <- scam(
  prob ~ s(loss, bs = "mpi", k = 15),
  data    = train_dat,
  weights = train_dat$n_at_risk,
  family  = binomial()
)


predict_escalation <- function(loss_vals, se = FALSE) {
  loss_clamped <- pmin(loss_vals, ME_threshold - 0.01)
  pred_link <- predict(gam_mod,
                       newdata = data.frame(loss = loss_clamped),
                       type = "link", se.fit = TRUE)
  prob  <- plogis(pred_link$fit)
  lower <- plogis(pred_link$fit - 1.96 * pred_link$se.fit)
  upper <- plogis(pred_link$fit + 1.96 * pred_link$se.fit)
  prob  <- pmin(pmax(prob,  0), 1)
  lower <- pmin(pmax(lower, 0), 1)
  upper <- pmin(pmax(upper, 0), 1)
  if (se) return(data.frame(prob = prob, lower = lower, upper = upper))
  return(prob)
}


# basal p (loss = 0%)
baseline_stats <- predict_escalation(0, se = TRUE)
sink('./NUMBERS/basal_p_escalation_model.txt')
cat(sprintf(
  "Baseline probability (loss = 0%%): %.2f%% (95%% CI: %.2f%% - %.2f%%)\n",
  baseline_stats$prob  * 100,
  baseline_stats$lower * 100,
  baseline_stats$upper * 100
))
sink()


###PLOT
ext_rates<-c(ext_rate_l,ext_rate_m,ext_rate_h)
plot_name<-c('low','interm','high')
for (plot_rep in 1:3){
  mu_best<-ext_rates[plot_rep]
  r_best  <- lambda - mu_best
  loss_seq <- seq(0, ME_threshold - 0.1, length.out = 400)
  
  esc_full_plot <- predict_escalation(loss_seq, se = TRUE) %>%
    mutate(
      loss       = loss_seq,
      delta_t    = log(1 - (loss / 100)) / r_best,
      year_actual = 2026 + delta_t,
      weight_convergence = pmax(0, (1 - prob)^2),
      lower_adj  = prob - (prob - lower) * weight_convergence,
      upper_adj  = pmin(1, prob + (upper - prob) * weight_convergence)
    )
  
  saturation_point <- esc_full_plot %>% filter(prob >= 0.999) %>% slice(1)
  end_year_dynamic <- if (nrow(saturation_point) > 0) {
    ceiling(saturation_point$year_actual / 100) * 100
  } else {
    10000
  }
  
  if (plot_rep==3){
    target_years <- c(2250, 2500, 2750, 3000)}
  else {
    target_years<-round(c(seq(2700,(end_year_dynamic*0.8),length.out=4)/10))*10
  }
  target_years <- target_years[target_years < end_year_dynamic]
  
  year_markers <- data.frame(year_actual = target_years) %>%
    mutate(
      delta_t = year_actual - 2026,
      loss    = (1 - exp(r_best * delta_t)) * 100,
      prob    = predict_escalation(loss)
    )
  
  legend_breaks <- unique(c(2026, 2500, end_year_dynamic))
  legend_breaks <- legend_breaks[legend_breaks <= end_year_dynamic]
  
  theme_common <- theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = 14),
      axis.title       = element_text(size = 16),
      legend.text      = element_text(size = 14),
      legend.title     = element_text(size = 14)
    )
  
  p_final <- ggplot(esc_full_plot %>% filter(year_actual <= end_year_dynamic),
                    aes(x = loss, y = prob)) +
    theme_bw() +
    theme_common +
    geom_ribbon(aes(ymin = lower_adj, ymax = upper_adj),
                fill = "grey90", alpha = 0.5) +
    geom_line(aes(color = year_actual), linewidth = 2.2, lineend = "round") +
    scale_color_viridis(
      option = "D",
      trans  = "log10",
      breaks = legend_breaks,
      labels = as.character(legend_breaks),
      limits = c(2026, end_year_dynamic),
      name   = "Year"
    ) +
    geom_point(data = year_markers, aes(x = loss, y = prob),
               size = 4, shape = 21, fill = "white", stroke = 1.2) +
    geom_text(data = year_markers,
              aes(x = loss, y = prob,
                  label = paste0(year_actual, "\n", round(prob * 100), "%")),
              vjust = -1.2, size = 4, fontface = "bold", lineheight = 0.8) +
    scale_y_continuous(
      name   = "Probability of escalation to mass extinction",
      labels = percent_format(),
      limits = c(0, 1.05),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      name   = "Diversity loss (%)",
      limits = c(0, max(esc_full_plot$loss[esc_full_plot$year_actual <= end_year_dynamic])),
      expand = c(0, 0)
    ) +
    theme(
      legend.key.height = unit(1.2, "cm"),
      plot.background   = element_rect(fill = "white", color = NA)
    )
  
  
  pdf(paste0("./FIGURES/escalation_final_",plot_name[plot_rep],".pdf"), width = 6, height = 5)
  print(p_final)
  dev.off()
}

###GAM diagnostics
summary(gam_mod)

# Residual Deviance (The error left in your model)
dev_res <- deviance(gam_mod)

# Null Deviance (The total error of a flat model)
# Since DevExplained = (Null - Res) / Null, then:
dev_null <- dev_res / (1 - 0.991) 

cat("Null Deviance:", dev_null, "\n")
cat("Residual Deviance:", dev_res, "\n")


