# [Trends in bird abundance differ among protected forests but not bird guilds](https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/eap.2377)

### Ecological Applications

### Jeffrey W. Doser, Aaron S. Weed, Elise F. Zipkin, Kathryn M. Miller, Andrew O. Finley

### Code/Data DOI: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4701477.svg)](https://doi.org/10.5281/zenodo.4701477)

### Please contact the first author for questions about the code or data: Jeffrey W. Doser (doserjef@msu.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:  
As birds continue to decline, monitoring population and community abundances is crucial for assessing effects on ecosystem functions and informing strategies to protect habitats and species. Improved monitoring and associated inferential tools to efficiently identify declining bird populations, particularly of rare or sparsely distributed species, is key to informed conservation and management across large spatio-temporal regions. We assess abundance trends for 106 bird species in a network of eight national park forests located within the northeast USA from 2006-2019 using a novel hierarchical model. We develop a multi-species, multi-region removal sampling model that shares information across species and parks to enable inference on rare species and sparsely sampled parks and to evaluate the effects of local forest structure. Trends in bird abundance over time varied widely across parks, but species showed similar trends within parks. Three parks (Acadia National Park and Marsh-Billings-Rockefeller and Morristown National Historic Parks (NHP)) decreased in bird abundance across all species, while three parks (Saratoga NHP and Roosevelt-Vanderbilt and Weir-Farm National Historic Sites) increased in abundance. Bird abundance peaked at medium levels of basal area and high levels of percent forest and forest regeneration, with percent forest having the largest effect. Variation in these effects across parks could be a result of differences in forest structural stage and diversity. Our novel hierarchical model enables uncertainty-quantified estimates of abundance at the network, park, guild, and species levels across a large spatio-temporal region. We found large variation in abundance trends across parks but not across bird guilds, suggesting that local forest condition may have a broad and consistent effect on the entire bird community within a given park. Management should target the three parks with overall decreasing trends in bird abundance to further identify what specific factors are driving observed declines across the bird community. Understanding how bird communities respond to local forest structure and other stressors (e.g., pest outbreaks, climate change) is crucial for informed and lasting management.

## [Published PDF](https://github.com/zipkinlab/Doser_etal_2021_EcoApps/blob/master/doser2021EcoApps.pdf)

## Code

1. [mr-ms-rs.R](https://github.com/doserjef/Doser_etal_2020_MR-MS-RS/blob/master/mr-ms-rs.R): Code to run hierarchical model in JAGS through R
2. [mr-ms-rs-jags.txt](https://github.com/doserjef/Doser_etal_2020_MR-MS-RS/blob/master/mr-ms-rs-jags.txt): JAGS code for multi-region, multi-species removal sampling (mr-ms-rs) model
3. [mr-ms-rs-summary.R](https://github.com/doserjef/Doser_etal_2020_MR-MS-RS/blob/master/mr-ms-rs-summary.R): Code to analyze the results of the mr-ms-rs model and code for figures in the manuscript

## Data

1. [mr-ms-rs-data.R](https://github.com/doserjef/Doser_etal_2020_MR-MS-RS/blob/master/mr-ms-rs-data.R): bird abundance data. Metadata is included in [mr-ms-rs.R](https://github.com/doserjef/Doser_etal_2020_MR-MS-RS/blob/master/mr-ms-rs.R)
