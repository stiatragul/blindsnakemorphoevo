# Alex Skeels
# 2023

# load packages
library(terra)# do everything with terra and sf instead of raster/sp/rgeos/rgdal
library(sf)
library(ape)
library(geiger)
library(phytools)
library(moments)
library(ggplot2)
library(abc)
library(ENMTools)
library(stringr)
library(apTreeshape)

# spatial
spp <- read.csv("./data/blindsnake_locality_for_AS.csv")

# read mcc tree
phy <- read.tree("./data/tr_anilios.tree")

# which are mismatched?
phy_not_spp <- phy$tip.label[which(!phy$tip.label %in% spp$species)]
spp_not_phy <- spp$species[which(!spp$species %in% phy$tip.label)]

# drop mismatched
spp <- spp[which(!spp$species %in% spp_not_phy ),]
phy <- drop.tip(phy, phy_not_spp)

# double check they match
all(spp$species %in% phy$tip.label)
all(phy$tip.label %in% spp$species)

# convert species occurences records into a raster
species_names <- unique(spp$species)
head(species_names)

# create a template raster which will set the resolution and extent of the rasters
raster_template <- rast(resolution=0.1, xmin=100, xmax=190, ymin=-50, ymax= 0)
values(raster_template) <- 0
crs(raster_template) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"# create some empty lists to fill up with rasters or spatial points
spp_stack <- list() # for spatial points

for(species_i in 1:length(species_names)){
  
  spp_i <- spp[which(spp$species==species_names[species_i]),]
  spp_i_sp <- SpatialPoints(spp_i[,c("longdec", "latdec")])
  spp_i <-  st_set_crs(st_as_sf(spp_i, coords=c("longdec", "latdec")),4326)
  spp_stack[[species_i]] <- spp_i
  if(species_i == 1){
    buffer_stack_01 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=1, fun=any)
    buffer_stack_03 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.3, fun=any)
    buffer_stack_05 <- rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.5, fun=any)
    mcp_stack       <-    rasterize(vect(adehabitatHR::mcp(spp_i_sp, percent = 90)), raster_template, field=1)  
  } else {
    buffer_stack_01 <- c(buffer_stack_01,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=1, fun=any))
    buffer_stack_03 <- c(buffer_stack_03,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.3, fun=any))
    buffer_stack_05 <- c(buffer_stack_05,rasterizeWin(x=data.frame(st_coordinates(spp_i), 1), y=raster_template, pars=0.5, fun=any))
    mcp_stack       <-    c(mcp_stack, rasterize(vect(adehabitatHR::mcp(spp_i_sp, percent = 90)), raster_template, field=1))  
  }
}


names(buffer_stack_01) <- species_names
names(buffer_stack_03) <- species_names
names(buffer_stack_05) <- species_names
names(mcp_stack)       <- species_names
names(spp_stack)       <- species_names


saveRDS(buffer_stack_01, file="output/Anilios_buffered_10.rds") 
saveRDS(buffer_stack_03, file="output/Anilios_buffered_03.rds")
saveRDS(buffer_stack_05, file="output/Anilios_buffered_05.rds") 
saveRDS(mcp_stack, file="output/Anilios_MCP.rds") 
saveRDS(spp_stack, file="output/Anilios_points.rds")   


scripts <-paste0("DREaD/", c("getSummaryStats_terrasf.R", "helperFunctions_terrasf.R","findSisters.R","summary_statsitics_functions_terrasf.R", "aoc.R")) 
sapply(scripts, source, echo=FALSE, verbose=FALSE)

#Get summary statistics of each clade for DREaD
emp_df <- data.frame(clade         = "Anilios",
                     ntips         = NA,
                     ROmean        = NA,
                     ROsd          = NA,
                     ROslope       = NA,
                     ROintercept   = NA,
                     ROskew        = NA,
                     ROkurtosis    = NA,
                     RO0           = NA,
                     RO50          = NA,
                     RO75          = NA,
                     RO90          = NA,
                     RO100         = NA,
                     TOmean        = NA,
                     TOsd          = NA,
                     asymmean      = NA,
                     asymsd        = NA,
                     asymslope     = NA,
                     asymintercept = NA,
                     RDmean        = NA,
                     RDsd          = NA,
                     RDintercept   = NA,
                     RDslope       = NA,
                     BIMOD50       = NA,
                     BIMODE75      = NA,
                     BIMOD90       = NA,
                     BIMOD100      = NA,
                     RSskew        = NA,
                     RSmean        = NA,
                     RSsd          = NA,
                     CI            = NA,
                     Beta          = NA,
                     Gamma         = NA,
                     SI            = NA,
                     ARCslope      = NA,
                     ARCint        = NA)


buffer_stack_01 <- readRDS(file="output/Anilios_buffered_10.rds") 
buffer_stack_03 <- readRDS(file="output/Anilios_buffered_03.rds")
buffer_stack_05 <- readRDS(file="output/Anilios_buffered_05.rds") 
mcp_stack <- readRDS(file="output/Anilios_MCP.rds") 


emp_df[1,] <- getSummaryStats(phy, "Anilios_buff01", buffer_stack_01)
emp_df[2,] <- getSummaryStats(phy, "Anilios_buff03", buffer_stack_03)
emp_df[3,] <- getSummaryStats(phy, "Anilios_buff05", buffer_stack_05)
emp_df[4,] <- getSummaryStats(phy, "Anilios_MCP", mcp_stack)


write.csv(emp_df, file=paste("output/Anilios_empirical_data.csv",sep=""), row.names=F)

## Simulated summary statistics for ABC

#Load and prepare DREaD simulation results
sim_df <- read.csv("DREaD/simulation_results_subset.csv",stringsAsFactors = F)

# get mode from which sims were generated 
speciation_modes <- as.character(sim_df$speciation_mode)
sim_df <- sim_df[,which(!colnames(sim_df) %in% "speciation_mode")]


#Model selection to infer the predominant geographic mode of speciation: ABC
emp_df_s <- emp_df[, which(colnames(emp_df) %in% colnames(sim_df))]
emp_df_s <- emp_df_s[,match(colnames(sim_df), colnames(emp_df_s))]

# get clade names
clades <- emp_df$clade

# create an empty list to put model results
abc_list <- list()

for(sub_clade_i in 1:length(clades)){
  
  # run ABC with logistic regression correction
  abc_log <- postpr(emp_df_s[sub_clade_i,], speciation_modes, sim_df, tol=.05, method="mnlogistic")
  
  # run ABC with lneural network correction
  abc_neu <- postpr(emp_df_s[sub_clade_i,], speciation_modes, sim_df, tol=.05, method="neuralnet", trace=F)
  
  # combine results into a list and store
  res <- cbind(abc_log$pred,abc_neu$pred)
  colnames(res) <- c("ABC_log", "ABC_nn")
  abc_list[[sub_clade_i]] <- round(res, 3)
  
}

# name list
names(abc_list) <- clades

# print results -> higher posterior probability means more likely predominant mode of speciation according to that algorithm (neural network or logistic regression)
print(abc_list)


saveRDS(abc_list, file="output/Anilios_ABC.rds")
abc_list <- readRDS( file="output/Anilios_ABC.rds")


abc_list <- lapply(abc_list, as.data.frame)
abc_list$Anilios_buff01$range="buff01"
abc_list$Anilios_buff03$range="buff03"
abc_list$Anilios_buff05$range="buff05"
abc_list$Anilios_MCP$range="MCP"
abc_list$Anilios_buff01$mode <- c("allopatric-vicariant", "allopatric-dispersal", "mixed", "parapatric", "sympatric")
abc_list$Anilios_buff03$mode <- c("allopatric-vicariant", "allopatric-dispersal", "mixed", "parapatric", "sympatric")
abc_list$Anilios_buff05$mode <- c("allopatric-vicariant", "allopatric-dispersal", "mixed", "parapatric", "sympatric")
abc_list$Anilios_MCP$mode <- c("allopatric-vicariant", "allopatric-dispersal", "mixed", "parapatric", "sympatric")

abc_df <- do.call(rbind, abc_list)

aggregate(abc_df$ABC_log, by=list(abc_df$mode), FUN=sum)[,2]+
aggregate(abc_df$ABC_nn, by=list(abc_df$mode), FUN=sum)[, 2]

cbPalette <- RColorBrewer::brewer.pal(5, "PuOr")


nn <- ggplot(abc_df, aes(x=range, y=ABC_nn, fill=mode)) + 
  labs(x="")+
  geom_bar(stat="identity") +
  ylab("Posterior (NN)") +
  theme_bw()+
  scale_fill_manual(values=cbPalette) +
  theme( panel.border = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         plot.background = element_rect(fill = 'white', colour = 'white'), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         panel.background = element_rect(fill = 'white', colour = 'white'))  

mnl <- ggplot(abc_df, aes(x=range, y=ABC_log, fill=mode)) + 
  labs(x="")+
  geom_bar(stat="identity") +
  ylab("Posterior (MnL)") +
  theme_bw()+
  scale_fill_manual(values=cbPalette) +
  theme( panel.border = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         plot.background = element_rect(fill = 'white', colour = 'white'), 
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         panel.background = element_rect(fill = 'white', colour = 'white'))  

library(gridExtra)
grid.arrange(nn, mnl)
