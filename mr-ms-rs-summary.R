# Author: Jeffrey W. Doser
# Code: analysis code to summarize the results of the multi-region, 
#       multi-species removal sampling model. 

# Empty memory
# rm(list = ls())


# Load packages -----------------------------------------------------------
library(coda)
library(jagsUI)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(patchwork)

# Load data ---------------------------------------------------------------
load("mr-ms-rs-data.R")
# Read in results ---------------------------------------------------------
# The results file 'mr-ms-rs-results.R' comes from the script 'mr-ms-rs.R'
# after running the model. Final results presented in Doser et al. (2021)
# are not included on GitHub as the file size is too large. Contact 
# Jeff Doser (doserjef@msu.edu) if file is desired. 
load("mr-ms-rs-results.R")
  # Parks: (1) ACAD (2) MABI (3) MIMA (4) MORR (5) ROVA 
         # (6) SAGA (7) SAGA (8) WEFA
# Number of iterations for posterior summary. 
n.iter <- 2500
# Save for use in derived quantities
mu.p.samples <- out$sims.list$mu.p[1:n.iter]
beta.1.p.samples <- out$sims.list$beta.1.p[1:n.iter]
beta.2.p.samples <- out$sims.list$beta.2.p[1:n.iter]
beta.3.p.samples <- out$sims.list$beta.3.p[1:n.iter]
beta.4.p.samples <- out$sims.list$beta.4.p[1:n.iter, , ]
mu.lambda.samples <- out$sims.list$mu.lambda[1:n.iter]
beta.1.lambda.samples <- out$sims.list$beta.1.lambda[1:n.iter]
beta.2.lambda.samples <- out$sims.list$beta.2.lambda[1:n.iter]
beta.3.lambda.samples <- out$sims.list$beta.3.lambda[1:n.iter]
beta.4.lambda.samples <- out$sims.list$beta.4.lambda[1:n.iter]
beta.5.lambda.samples <- out$sims.list$beta.5.lambda[1:n.iter]
mu.p.r.samples <- out$sims.list$mu.p.r[1:n.iter, ]
beta.1.p.r.samples <- out$sims.list$beta.1.p.r[1:n.iter, ]
beta.2.p.r.samples <- out$sims.list$beta.2.p.r[1:n.iter, ]
beta.3.p.r.samples <- out$sims.list$beta.3.p.r[1:n.iter, ]
beta.4.p.r.samples <- out$sims.list$beta.4.p.r[1:n.iter, ]
mu.lambda.r.samples <- out$sims.list$mu.lambda.r[1:n.iter, ]
beta.1.lambda.r.samples <- out$sims.list$beta.1.lambda.r[1:n.iter, ]
beta.2.lambda.r.samples <- out$sims.list$beta.2.lambda.r[1:n.iter, ]
beta.3.lambda.r.samples <- out$sims.list$beta.3.lambda.r[1:n.iter, ]
beta.4.lambda.r.samples <- out$sims.list$beta.4.lambda.r[1:n.iter, ]
beta.5.lambda.r.samples <- out$sims.list$beta.5.lambda.r[1:n.iter, ]
mu.p.sp.samples <- out$sims.list$mu.p.sp[1:n.iter, , ]
beta.1.p.sp.samples <- out$sims.list$beta.1.p.sp[1:n.iter, , ]
beta.2.p.sp.samples <- out$sims.list$beta.2.p.sp[1:n.iter, , ]
beta.3.p.sp.samples <- out$sims.list$beta.3.p.sp[1:n.iter, , ]
beta.4.p.sp.samples <- out$sims.list$beta.4.p.sp[1:n.iter, , ]
mu.lambda.sp.samples <- out$sims.list$mu.lambda.sp[1:n.iter, , ]
beta.1.lambda.sp.samples <- out$sims.list$beta.1.lambda.sp[1:n.iter, , ]
beta.2.lambda.sp.samples <- out$sims.list$beta.2.lambda.sp[1:n.iter, , ]
beta.3.lambda.sp.samples <- out$sims.list$beta.3.lambda.sp[1:n.iter, , ]
beta.4.lambda.sp.samples <- out$sims.list$beta.4.lambda.sp[1:n.iter, , ]
beta.5.lambda.sp.samples <- out$sims.list$beta.5.lambda.sp[1:n.iter, , ]
bp.1.samples <- out$sims.list$bp.1[1:n.iter]
N.samples <- out$sims.list$N[1:n.iter, , , , ]


# Network-level Trend -----------------------------------------------------
summary(mcmc(beta.1.lambda.samples))

# Park-level Trends -------------------------------------------------------
summary(mcmc(beta.1.lambda.r.samples))

# Species-level trend summary ---------------------------------------------

# Get trends for each species at all parks where species is observed
my.years <- sort(unique(c(year.regions)))
sp.trends <- list()
sp.year.estimates <- list()
regions <- c("ACAD", "MABI", "MIMA", "MORR", "ROVA", "SAGA", "SARA", "WEFA")
all.years <- 2006:2019
for (j in 1:n) {
  print(paste("Species: ", j, " out of ", n, sep = ''))
  curr.sp <- apply(r.sp, 1, function(a) {which(a == sp[j])})
  names(curr.sp) <- regions
  curr.sp <- unlist(curr.sp)
  names(curr.sp) <- substr(names(curr.sp), 1, 4)
  region.indx <- which(regions %in% names(curr.sp))
  my.regions <- regions[region.indx]
  my.year.reg <- year.regions[region.indx]
  curr.R <- length(region.indx)

  # Get posterior samples of regression coefficients for species of interest
  # The plus 1 is ito include the mean for the NETN
  curr.beta.1.samples <- matrix(NA, n.iter, curr.R + 1)
  for (i in 1:length(region.indx)) {
    curr.year.ind <- which(all.years %in% year.regions[region.indx[i], ])
    curr.beta.1.samples[, i] <- beta.1.lambda.sp.samples[, curr.sp[i], region.indx[i]]
  }
  beta.1.netn.samples <- apply(matrix(curr.beta.1.samples[, -(curr.R + 1)], n.iter, curr.R),
			       1, mean, na.rm = TRUE)
  curr.beta.1.samples[, curr.R + 1] <- beta.1.netn.samples

  # Get year effect estimate for each species
  year.meds <- apply(curr.beta.1.samples, 2, median, na.rm = TRUE)
  names(year.meds) <- c(my.regions, 'NETN')
  year.low <- apply(curr.beta.1.samples, 2, quantile, na.rm = TRUE, probs = c(0.025))
  year.high <- apply(curr.beta.1.samples, 2, quantile, na.rm = TRUE, probs = c(0.975))
  names(year.high) <- c(my.regions, 'NETN')
  names(year.low) <- c(my.regions, 'NETN')
  year.vals <- rbind(year.low, year.meds, year.high)
  tmp <- data.frame(year.vals)
  tmp$type <- c('lowCI', 'median', 'highCI')
  sp.year.estimates[[j]] <- tmp
}

sp.years.df <- bind_rows(sp.year.estimates, .id = 'column_label')
sp.years.df$species <- rep(sp, each = 3)
sp.years.df <- sp.years.df %>% select(species, type, ACAD, MABI, MIMA, MORR, ROVA, SARA, SAGA, WEFA, NETN)

sp.years.long <- pivot_longer(sp.years.df, ACAD:NETN, names_to = 'Park', values_to = 'val')
sigs <- matrix(NA, n, R + 1)

# Determine number of significant trends
parks <- c(regions, 'NETN')
for (i in 1:n) {
  for (r in 1:(R+1)) {
    curr.dat <- sp.years.long %>%
      filter(species == sp[i], Park == parks[r])
    if (sum(is.na(curr.dat$val)) == 0) {
      low <- curr.dat %>% filter(type == 'lowCI') %>% pull(val)
      high <- curr.dat %>% filter(type == 'highCI') %>% pull(val)
      sigs[i, r] <- ifelse((0 > low) & (0 < high), 0, 1)
    }    
  }  
}

# For plots in Figure 3
plot.dat <- sp.years.long %>%
  filter(type == 'median') %>%
  select(species, Park, val) %>%
  mutate(significant = c(t(sigs)),
	 Park = factor(Park, levels = c('ACAD', 'MABI', 'MIMA', 'MORR', 'ROVA',
					'SAGA', 'SARA', 'WEFA', 'NETN')),
	 species = factor(species, levels = rev(levels(species))))
plot.dat$missing <- ifelse(is.na(plot.dat$val), '-', '')
plot.dat <- plot.dat %>%
  mutate(group = factor(ifelse((val < 0) & (significant == 1), '-', 
			       ifelse((val > 0) & (significant == 1), '+', 
				      '0')), levels = c('-', '0', '+')), 
	 sign = factor(ifelse(val < 0, '-', '+'), levels = c('-', '+')))

# Example plot for bar plot in Figure 3
cols <- c('#ff8f80', '#c1e4f7')
curr.park <- 'ACAD'
sig.vals <- plot.dat %>%
  filter(Park == curr.park, !is.na(val)) %>%
  group_by(group, .drop = FALSE) %>%
  summarize(n.vals = n())
plot.lims <- plot.dat %>%
  filter(Park == curr.park, !is.na(val)) %>%
  group_by(sign) %>%
  summarize(n.vals = n())
plot.dat %>% 
  filter(Park == curr.park, !is.na(val)) %>%
  ggplot(aes(x = sign, fill = sign)) + 
    geom_bar(position = 'dodge', show.legend = FALSE, col = 'black', size = 2) + 
    scale_fill_manual(values = cols) + 
    theme_bw() + 
    labs(x = 'Trend', y = '# of Species', title = 'Acadia National Park (ACAD)') +
    theme(text = element_blank(), 
	  axis.ticks.y = element_blank(), 
	  axis.ticks.x = element_blank()) + 
    ylim(0, 100) 
# Other plots are analagous, not included here to reduce clutter. 

# Species Richness --------------------------------------------------------

rich.samples <- array(NA, c(n.iter, R, n.years))

for (a in 1:n.iter) {
  print(a)
  for (r in 1:R) {
    for (t in 1:year.by.regions[r]) {
       curr.vals <- ifelse(N.samples[a, r, , 1:J[r], t] > 0, 1,
			   N.samples[a, r, , 1:J[r], t])
       tmp.1 <- apply(curr.vals, 1, function(a) {ifelse(sum(a, na.rm = TRUE) > 0, 1, 0)})
       rich.samples[a, r, t] <- sum(tmp.1, na.rm = TRUE)
    }
  }
}

rich.ordered.samples <- array(NA, c(n.iter, R + 1, n.years))
all.years <- 2006:2019
for (r in 1:R) {
  curr.years <- which(all.years %in% year.regions[r, ])
  rich.ordered.samples[, r, curr.years] <- rich.samples[, r, 1:year.by.regions[r]]
}
rich.netn.samples <- apply(rich.ordered.samples[, -(R + 1), ], 
			   c(1, 3), mean, na.rm = TRUE)
rich.ordered.samples[, R + 1, ] <- rich.netn.samples
rich.med <- as.vector(t(apply(rich.ordered.samples, c(2, 3), mean, na.rm = TRUE)))
rich.lowest <- as.vector(t(apply(rich.ordered.samples, c(2, 3), 
			      quantile, probs = 0.025, na.rm = TRUE)))
rich.low <- as.vector(t(apply(rich.ordered.samples, c(2, 3), 
			   quantile, probs = 0.25, na.rm = TRUE)))
rich.high <- as.vector(t(apply(rich.ordered.samples, c(2, 3), 
			    quantile, probs = 0.75, na.rm = TRUE)))
rich.highest <- as.vector(t(apply(rich.ordered.samples, c(2, 3), 
			       quantile, probs = 0.975, na.rm = TRUE)))
rich.med <- rich.med[(rich.med != 0) & (!is.na(rich.med))]
rich.low <- rich.low[!is.na(rich.low)]
rich.lowest <- rich.lowest[!is.na(rich.lowest)]
rich.high <- rich.high[!is.na(rich.high)]
rich.highest <- rich.highest[!is.na(rich.highest)]
years.tracked <- as.vector(t(year.regions))
years.tracked <- years.tracked[!is.na(years.tracked)]
# NETN is observed every year
years.tracked <- c(years.tracked, sort(unique(years.tracked))) 
rich.plot.dat <- data.frame(years = years.tracked,
                           region = factor(rep(c(regions, 'NETN'), 
					c(year.by.regions, n.years)),
					   levels = c(regions, 'NETN')),
                           lowest = rich.lowest,
                           low = rich.low,
                           richness = rich.med,
                           high = rich.high,
                           highest = rich.highest)


pdf("acadRichness.pdf", width = 3, height = 3)
rich.plot.dat %>%
  filter(region == 'ACAD') %>%
ggplot(aes(x = years, y = richness)) +
geom_point(size = 5) +
geom_line(size = 2) +
theme_void(base_size = 40) +
geom_ribbon(aes(x = years, ymin = lowest, ymax = highest),
            fill = "gray80", alpha = 0.4, col = NA) +
theme(axis.ticks.x = element_blank(), 
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()) +  
scale_y_continuous(limits = c(10, 70)) 
dev.off()

# Guild-level Trends ------------------------------------------------------
canopy.nest <- guilds %>%
  filter(Response_Guild == 'CanopyNester') %>%
  pull(AOU_Code)
ground.nest <- guilds %>%
  filter(Response_Guild == 'ForestGroundNester') %>%
  pull(AOU_Code)
shrub.nest <- guilds %>%
  filter(Response_Guild == 'ShrubNester') %>%
  pull(AOU_Code)
forest.obligate <- guilds %>%
  filter(Response_Guild == 'InteriorForestObligate') %>%
  pull(AOU_Code)
open.ground.nest <- guilds %>%
  filter(Response_Guild == 'OpenGroundNester') %>%
  pull(AOU_Code)
forest.gen <- guilds %>%
  filter(Response_Guild == 'ForestGeneralist') %>%
  pull(AOU_Code)
omnivore <- guilds %>% 
  filter(Response_Guild == 'Omnivore') %>%
  pull(AOU_Code)
bark.prober <- guilds %>% 
  filter(Response_Guild == 'BarkProber') %>%
  pull(AOU_Code)
ground.gleaner <- guilds %>% 
  filter(Response_Guild == 'GroundGleaner') %>%
  pull(AOU_Code)
upper.can.for <- guilds %>% 
  filter(Response_Guild == 'UpperCanopyForager') %>%
  pull(AOU_Code)
lower.can.for <- guilds %>%
  filter(Response_Guild == 'LowerCanopyForager') %>%
  pull(AOU_Code)
exotic <- guilds %>%
  filter(Response_Guild == 'Exotic') %>%
  pull(AOU_Code)
resident <- guilds %>%
  filter(Response_Guild == 'Resident') %>%
  pull(AOU_Code)
single.brooded <- guilds %>%
  filter(Response_Guild == 'SingleBrooded') %>%
  pull(AOU_Code)
nest.pred <- guilds %>% 
  filter(Response_Guild == 'NestPredator_BroodParasite') %>%
  pull(AOU_Code)
temp.mig <- guilds %>% 
  filter(Response_Guild == 'TemperateMigrant') %>%
  pull(AOU_Code)

guild.list <- list(omnivore, bark.prober, ground.gleaner, upper.can.for,
		   lower.can.for, exotic, resident, single.brooded, nest.pred,
		   temp.mig, canopy.nest, shrub.nest, ground.nest,
		   forest.obligate, forest.gen, open.ground.nest)

n.guilds <- 16
beta.1.guilds.samples <- array(NA, dim = c(n.iter, R, n.guilds))
beta.2.guilds.samples <- array(NA, dim = c(n.iter, R, n.guilds))
beta.3.guilds.samples <- array(NA, dim = c(n.iter, R, n.guilds))
beta.4.guilds.samples <- array(NA, dim = c(n.iter, R, n.guilds))
beta.5.guilds.samples <- array(NA, dim = c(n.iter, R, n.guilds))

for (i in 1:n.guilds) {
  print(i)
  curr.guilds <- as.character(guild.list[[i]])
  for (r in 1:R) {
    curr.reg <- regions[r]
    reg.sp <- as.matrix(r.sp[r, ])
    curr.sp <- which(reg.sp %in% curr.guilds)
    curr.beta.1.sp.samples <- matrix(beta.1.lambda.sp.samples[, curr.sp, r], nrow = n.iter,
                                     ncol = length(curr.sp))
    beta.1.guilds.samples[, r, i] <- apply(curr.beta.1.sp.samples, 1, mean)
    beta.2.guilds.samples[, r, i] <- apply(matrix(beta.2.lambda.sp.samples[, curr.sp, r],
                                                  nrow = n.iter), 1, mean)
    beta.3.guilds.samples[, r, i] <- apply(matrix(beta.3.lambda.sp.samples[, curr.sp, r],
                                                  nrow = n.iter), 1, mean)
    beta.4.guilds.samples[, r, i] <- apply(matrix(beta.4.lambda.sp.samples[, curr.sp, r],
                                                  nrow = n.iter), 1, mean)
    beta.5.guilds.samples[, r, i] <- apply(matrix(beta.5.lambda.sp.samples[, curr.sp, r],
                                                  nrow = n.iter), 1, mean)
  }
}


# For table of guild effects at the network level
beta.1.guild.netn <- apply(beta.1.guilds.samples, c(1, 3), mean, na.rm = TRUE)
beta.2.guild.netn <- apply(beta.2.guilds.samples, c(1, 3), mean, na.rm = TRUE)
beta.3.guild.netn <- apply(beta.3.guilds.samples, c(1, 3), mean, na.rm = TRUE)
beta.4.guild.netn <- apply(beta.4.guilds.samples, c(1, 3), mean, na.rm = TRUE)
beta.5.guild.netn <- apply(beta.5.guilds.samples, c(1, 3), mean, na.rm = TRUE)

summary(mcmc(beta.4.guild.netn))$quantiles
summary(mcmc(beta.5.guild.netn))$quantiles
summary(mcmc(beta.3.guild.netn))$quantiles
summary(mcmc(beta.2.guild.netn))$quantiles
summary(mcmc(beta.1.guild.netn))$quantiles

# Create Figure 5 ---------------------

guilds.ordered <- c('Omnivore', 'BarkProber', 'GroundGleaner', 'UpperCanopyForager', 
		    'LowerCanopyForager', 'Exotic', 'Resident', 'SingleBrooded', 
		    'NestPredator_BroodParasite', 'TemperateMigrant', 
		    'CanopyNester', 'ShrubNester', 'ForestGroundNester', 
		    'InteriorForestObligate', 'ForestGeneralist', 'OpenGroundNester')
guilds.region <- list()

for (r in 1:R) {
  curr.sp <- unlist(lapply(r.sp[r, ], as.character))
  curr.guilds <- guilds %>%
    filter(AOU_Code %in% curr.sp) %>%
    pull(Response_Guild) %>%
    unique()
  guilds.region[[r]] <- which(guilds.ordered %in% curr.guilds)    
}


year.guild.meds <- matrix(NA, R, n.guilds)
year.guild.low <- matrix(NA, R, n.guilds)
year.guild.high <- matrix(NA, R, n.guilds)
for (r in 1:R) {
  year.guild.meds[r, guilds.region[[r]]] <- 
	  summary(mcmc(beta.1.guilds.samples[, r, guilds.region[[r]]]))$quantiles[, 3]
  year.guild.low[r, guilds.region[[r]]] <-
	summary(mcmc(beta.1.guilds.samples[, r, guilds.region[[r]]]))$quantiles[, 1]
  year.guild.high[r, guilds.region[[r]]] <- 
	  summary(mcmc(beta.1.guilds.samples[, r, guilds.region[[r]]]))$quantiles[, 5]
}

year.meds.df <- data.frame(year.guild.meds)
names(year.meds.df) <- c('Omnivore', 'Bark Prober', 'Ground Gleaner', 'Upper Canopy Forager', 
			 'Lower Canopy Forager', 'Exotic', 'Resident', 'Single Brooded', 
			 'Nest Predator/Brood Parasite', 'Temperate Migrant', 
			 'Canopy Nester', 'Shrub Nester', 'Forest Ground Nester', 
			 'Interior Forest Obligate', 'Forest Generalist', 'Open Ground Nester')
year.meds.df$park <- factor(regions, levels = c('MABI', 'ACAD', 'MORR', 'SAGA', 'MIMA', 
						'ROVA', 'SARA', 'WEFA'))
dat.long <- pivot_longer(year.meds.df, 1:16, 
			 names_to = 'Guild', values_to = 'val')

year.low.df <- data.frame(year.guild.low)
names(year.low.df) <- c('Omnivore', 'Bark Prober', 'Ground Gleaner', 'Upper Canopy Forager', 
			 'Lower Canopy Forager', 'Exotic', 'Resident', 'Single Brooded', 
			 'Nest Predator/Brood Parasite', 'Temperate Migrant', 
			 'Canopy Nester', 'Shrub Nester', 'Forest Ground Nester', 
			 'Interior Forest Obligate', 'Forest Generalist', 'Open Ground Nester')
year.low.df$park <- factor(regions, levels = c('MABI', 'ACAD', 'MORR', 'SAGA', 'MIMA', 
						'ROVA', 'SARA', 'WEFA'))
low.long <- pivot_longer(year.low.df, 1:16, 
			 names_to = 'Guild', values_to = 'low')

year.high.df <- data.frame(year.guild.high)
names(year.high.df) <- c('Omnivore', 'Bark Prober', 'Ground Gleaner', 'Upper Canopy Forager', 
			 'Lower Canopy Forager', 'Exotic', 'Resident', 'Single Brooded', 
			 'Nest Predator/Brood Parasite', 'Temperate Migrant', 
			 'Canopy Nester', 'Shrub Nester', 'Forest Ground Nester', 
			 'Interior Forest Obligate', 'Forest Generalist', 'Open Ground Nester')
year.high.df$park <- factor(regions, levels = c('MABI', 'ACAD', 'MORR', 'SAGA', 'MIMA', 
						'ROVA', 'SARA', 'WEFA'))
high.long <- pivot_longer(year.high.df, 1:16, 
			 names_to = 'Guild', values_to = 'high')
dat.long$low <- low.long$low
dat.long$high <- high.long$high
dat.long <- dat.long %>%
  mutate(sig = ifelse(((low < 0) & (high < 0)) | ((low > 0) & (high > 0)), '*', ''))

png("guildSummarySig.png", width = 7, height = 7, res = 100, units = 'in')
ggplot(dat.long, aes(x = park, y = Guild, fill = val)) + 
  geom_raster() +
  scale_fill_gradient2(midpoint = 0, low = 'red', mid = 'white', high = 'blue', 
		       na.value = 'gray') + 
  geom_text(aes(label = sig)) + 
  theme_bw() + 
  labs(x = 'Park', y = 'Guild', fill = 'Trend') + 
  theme(axis.ticks.y = element_blank(), 
	axis.ticks.x = element_blank())
dev.off()

# Covariate Effects -------------------------------------------------------

# Regeneration
summary(mcmc(beta.2.lambda.r.samples))
# Percent Forest
summary(mcmc(beta.3.lambda.r.samples))
# Basal Area
summary(mcmc(beta.4.lambda.r.samples))
# Basal Area^2
summary(mcmc(beta.5.lambda.r.samples))

# Detection Effects -------------------------------------------------------

# Day
summary(mcmc(beta.1.p.r.samples))
# Day^2
summary(mcmc(beta.2.p.r.samples))
# Time
summary(mcmc(beta.3.p.r.samples))
