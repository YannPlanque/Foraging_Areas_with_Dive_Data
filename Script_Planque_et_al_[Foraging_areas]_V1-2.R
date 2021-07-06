##--------------------------------------------------------------------------------------------------------
## SCRIPT : Identification of seal foraging areas from likely foraging dives selected with a vertical approach
## Specific content : - Dive data analysis:
##                       => Calculate individual Minimum Cost of Transport Speed (MCTS)
##                       => Determine dive shape with the Time Allocation at Depth (TAD) index (Fedak et al., 2001).
##                       => Select likely foraging dives using two dive criteria:
##                          dive shape (TAD) and vertical descent speed
##                          (vertical approach by Planque et al., 2020)
##                    - Spatialise foraging areas with likely foraging dives (here with Kernel contours)
##
## As part of : 
##    Planque Y, Spitz J, Authier M, Guillou G, Vincent C, Caurant F.
##    Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus)
##     at the southern limit of their European range (Eastern English Channel).
##     Ecology and Evolution. 2021;00:1– 22. https://doi.org/10.1002/ece3.7739
##
## Using the vertical approach presented by :
##    Planque, Y., Huon, M., Caurant, F., Pinaud, D., Vincent, C., 2020.
##     Comparing the horizontal and vertical approaches used to identify foraging areas
##     of two diving marine predators. Marine Biology 167, 25.
##     https://doi.org/10.1007/s00227-019-3636-8
##
##
## Author : Yann Planque(1)*
## Affiliation : 
##    (1) Centre d'Etudes Biologiques de Chize (CEBC, UMR 7372 CNRS - La Rochelle Universite), La Rochelle, France
##
## Contact* : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
##
## First publication on GitHub : 2021-04-14
## Last update : 2021-07-06 (Version 1.2)
##
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##---------------------------------------------------------------------------------------------------------


### 0 // Packages ##########################################################################################

lapply(c("dplyr", "tidyr", "ggplot2", "ggpubr", "gridExtra",      # for plots
         "KernSmooth",                                            # for Kernels
         "rgdal", "sp", "sf", "grDevices", "maptools", "ggsn"),   # for spatial
       library, character.only=TRUE)

###########################################################################################################


### 0 // Useful functions #################################################################################

  ## 1 / Dive shape 

### Function "MCTS_ind" to calculate the MCTS with dive data
# We recommand to calculate MCTS at the individual level to consider interindividual differences.
# Based on the calculation of TAD index describing dive shape (Fedak et al., 2001; see the function TAD_index below).
# Method: presented by Vincent et al. (2016) (https://doi.org/10.1093/icesjms/fsw102):
#         "The best S value was selected for each seal when the second derivative of the
#         number of 0.5 < TAD < 1.0 came to zero (or changed from a negative to a positive value).
# Parameters: 3 dive parameters needed (colnames in the data can be set with the function):
#                 - dive duration (default "DIVE_DUR"), in s
#                 - maximum depth (default "MAX_DEP"), in m
#                 - percent area (default "PERCENT_AREA"), in %
#             Dive threshold = the limit of dive detection (default is 1.5 m depth)
#             Test MCTS from min (default: 1.0) to max (default: 3.0), by a specific value (default: 0.1)
# Packages needed: ggplot2, ggpubr, gridExtra.
# Data input: dive data with needed parameters.
# Output: the estimated MCTS (in m.s-1).

MCTS_ind <- function(data, DIVE_DUR = "DIVE_DUR", MAX_DEP = "MAX_DEP", PERCENT_AREA = "PERCENT_AREA", 
                     MCTS_min = 1.0, MCTS_max = 3.0, MCTS_by = 0.1, Dive_threshold = 1.5,
                     return_plot_obj = F, show_plot = T){
  
  S <- seq(MCTS_min, MCTS_max, MCTS_by) # S = Vertical travel speeds (i.e. MCTS) tested in TAD calculation
  
  # 1// Calculate the number of dives with TAD > 0.5 & TAD < 1.0 for each MCTS tested (noted "S" value below)
  MCTS_selec <- NULL
  for (i in 1:length(S)) {  
    
    S_test <- S[i]
    
    # TAD index (see the formula in Vincent et al. 2016)
    TAD <- (((data[[PERCENT_AREA]]/100) * (data[[DIVE_DUR]]) * (data[[MAX_DEP]]-Dive_threshold)) -
              (((data[[MAX_DEP]] - Dive_threshold) ^ 2) / S_test)) /
      ((data[[DIVE_DUR]]) * (data[[MAX_DEP]] - Dive_threshold)-
         2 * (((data[[MAX_DEP]] - Dive_threshold)^2) / S_test))
    
    nb_TAD <- sum(TAD > 0.5 & TAD < 1.0, na.rm = T)   # Number of dives with TAD conditions
    
    MCTS_ind <- data.frame(S=S_test, NTAD05=nb_TAD)   # summarise informations in MCTS_selec for each S-MCTS tested 
    MCTS_selec <- rbind(MCTS_selec,MCTS_ind)   
  }
  
  
  # 1// Identify the best S-MCTS tested (method in Vincent et al., 2016)
  # A// Calculate first and second derivatives between successive points with diff()
  MCTS_selec$First_deriv <- c(diff(MCTS_selec$NTAD05),0)  
  MCTS_selec$Second_deriv <- c(diff(diff(MCTS_selec$NTAD05)),0,0)  
  
  # B// Select the best S-MCTS: When the number of NTAD05 according to the S value reaches an asymptote (on plot). 
  #     i.e. when when the second derivative of the number of dives
  #     with 0.5 < TAD < 1.0 came to zero (or changed from a negative to a positive value).
  # 
  # Note : when the result give a sigmoid, the second derivative will be positive > negative > positive. 
  #        In this case, the selected S is the first S with second derivative becoming positive (i.e. negative > positive),
  #        and not the first positive value
  
  First_negativ_sec_der <- MCTS_selec[MCTS_selec$Second_deriv < 0,][1,1]  # First negative val
  
  # Select the first S with the first Second_deriv to become positive (from a previous negative value)
  MCTS_S_selected <- MCTS_selec[MCTS_selec$Second_deriv >= 0 &
                                  MCTS_selec$S > First_negativ_sec_der,][1,1] 
  
  
  ## PREPARATION OF FIGURE - second derivative
  # Data in table
  table_summary <- MCTS_selec[,c(1,2,4)]  
  names(table_summary) <- c("S","Nb_dives","Sec_deriv")  
  
  table_summary$MCTS <- ifelse(table_summary$S == MCTS_S_selected, 1, 0)
  
  range_N <- max(table_summary$Nb_dives) - min(table_summary$Nb_dives)
  Min_N <- min(table_summary$Nb_dives)
  range_sec <- max(table_summary$Sec_deriv) - min(table_summary$Sec_deriv)
  Min_sec <- min(table_summary$Sec_deriv)
  
  Plot_NB <-
    table_summary %>%
    ggplot(aes(S, Nb_dives)) +
    geom_segment(data = table_summary[table_summary$MCTS == 1,],
                 aes(x = MCTS_S_selected, xend = MCTS_S_selected,
                     y = Min_N - 0.065*range_N, yend = Nb_dives),
                 linetype = "dashed", size=1.4, color="chartreuse3") +
    geom_line(color="grey", size=1.15) +
    geom_point(shape=20, size=4.5) +
    geom_text(aes(label = S), vjust = -1, nudge_y = 0.5, size=4.5,
              color="royalblue3", fontface="bold") +
    geom_text(aes(x = MCTS_S_selected+0.2, y = Min_N,
                  label = paste0("MCTS = ", MCTS_S_selected)),
              size=5.5, color="chartreuse3") +
    labs(x="S", y="N dives with 0.5 < TAD < 1") +
    ylim(min(table_summary$Nb_dives) - 0.065*range_N, 
         max(table_summary$Nb_dives) + 0.08*range_N) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(size = 12, color="black"),
          axis.text.y = element_text(size = 12, color="black", angle=90, hjust = 0.5),
          axis.title = element_text(size = 14, color="black", face = "bold")) 
  
  Plot_sec_deriv <-
    table_summary %>%
    ggplot(aes(S, Sec_deriv)) +
    geom_hline(yintercept = 0, size=1, color="red2") +
    geom_segment(data = table_summary[table_summary$MCTS == 1,],
                 aes(x = MCTS_S_selected, xend = MCTS_S_selected,
                     y = Min_sec - 0.065*range_sec, yend = Sec_deriv),
                 linetype = "dashed", size=1.4, color="chartreuse3") +
    geom_line(color="grey", size=1.15) +
    geom_point(shape=20, size=4.5) +
    geom_text(aes(label = S), vjust = -1, nudge_y = 0.5, size=4.5,
              color="royalblue3", fontface="bold") +
    geom_text(aes(x = MCTS_S_selected+0.2, y = Min_sec,
                  label = paste0("MCTS = ", MCTS_S_selected)),
              size=5.5, color="chartreuse3") +
    labs(x="S", y="Second derivative") +
    ylim(min(table_summary$Sec_deriv) - 0.065*range_sec, 
         max(table_summary$Sec_deriv) + 0.08*range_sec) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(size = 12, color="black"),
          axis.text.y = element_text(size = 12, color="black", angle=90, hjust = 0.5),
          axis.title = element_text(size = 14, color="black", face = "bold")) 
  
  Plot_table <-
    ggplot(data.frame()) + geom_point() +
    xlim(-1, 1) + ylim(-1, 1) +
    annotation_custom(gridExtra::tableGrob(table_summary, rows=NULL, 
                                           theme=gridExtra::ttheme_default(base_size = 9))) +
    theme_void()
  
  Final_plot <-
    ggpubr::ggarrange(ggpubr::ggarrange(Plot_NB, Plot_sec_deriv, align = "v", ncol=1),
                      Plot_table, ncol = 2, widths = c(3,1))
  
  if(isTRUE(return_plot_obj)){
    if(isTRUE(show_plot)){
      print(Final_plot)}else{}
    return(Final_plot)
  }else{
    if(isTRUE(show_plot)){
      print(Final_plot)}else{}
    return(MCTS_S_selected)
  }
  
}


### Function "TAD_index" to calculate the Time Allocation at Depth (TAD) index with a specific MCTS.
# TAD index describes dive shape. TAD values mostly vary between 0 and 1.
# 0 = most of diving time close to the surface. 0.5 = V-shaped dive. 1 = U-shaped dive. TAD > 1 for U-shaped dives with a speed > the MCTS.
# TAD index developped by Fedak et al. (2001) (https://doi.org/10.1111/j.1748-7692.2001.tb00982.x)
# Used formula in Vincent et al. (2016) (https://doi.org/10.1093/icesjms/fsw102)
# Parameters: 4 are needed in data columns: DIVE_DUR (in s), MAX_DEP (in m), PERCENT_AREA (in %), MCTS (in m.s-1)
#             (indicate the column names if different from the default ones)
#             Dive threshold = the limit of dive detection (default is 1.5 m depth)
# Data input: dive data with needed parameters.
# Output: a TAD value for each dive (i.e. one value in the order of rows in input).

TAD_index <- function(data, DIVE_DUR = "DIVE_DUR", MAX_DEP = "MAX_DEP",
                      PERCENT_AREA = "PERCENT_AREA", MCTS = "MCTS", Dive_threshold = 1.5){
  
  if(is.numeric(MCTS)){   # If the MCTS is direcly specified
    Tad_values <- (((data[[PERCENT_AREA]]/100) * (data[[DIVE_DUR]]) * (data[[MAX_DEP]] - Dive_threshold)) -
                     (((data[[MAX_DEP]] - Dive_threshold) ^ 2) / MCTS)) /
      ((data[[DIVE_DUR]]) * (data[[MAX_DEP]] - Dive_threshold) -
         2 * (((data[[MAX_DEP]] - Dive_threshold) ^ 2) / MCTS))
    
  }else{    # If the MCTS is in column 
    Tad_values <- (((data[[PERCENT_AREA]]/100) * (data[[DIVE_DUR]]) * (data[[MAX_DEP]] - Dive_threshold)) -
                     (((data[[MAX_DEP]] - Dive_threshold) ^ 2) / data[[MCTS]])) /
      ((data[[DIVE_DUR]]) * (data[[MAX_DEP]] - Dive_threshold) -
         2 * (((data[[MAX_DEP]] - Dive_threshold) ^ 2) / data[[MCTS]]))
    
  }
  
  # Zeros if the animal spent all of this time at the depth of 'Dive_threshold'. In this case, TAD = 0.
  Tad_values <- ifelse(is.infinite(Tad_values), 0, Tad_values)
  Tad_values <- ifelse(is.na(Tad_values), 0, Tad_values)
  Tad_values <- ifelse(is.nan(Tad_values), 0, Tad_values)

  return(Tad_values)
  
}


  ## 2 / Kernel contours 

### Function calcHR() proposed by Fieberg (2014) (http://dx.doi.org/10.13020/D6G59W)
# Function to calculate Home Range, and needed to calculate Kernel contours afterwards (cf. function "Kernel_polyg_fast")

calcHR <- function(est, p = c(0.5, 0.95)){
  
  #  p = quantiles for estimation and display purposes
  
  x<-est$x1
  y<-est$x2
  z<-est$fhat
  
  # Calculate HR size
  #  Get size of grid
  dx<-x[2]-x[1]
  dy<-y[2]-y[1]
  
  # Estimate total volume and make sure =1 or close to
  totp<-sum(z)*dx*dy
  
  if(round(totp,4)<0.95){
    print("ERROR, total probability not equal to 1")
    print(paste("totp = ",totp))}
  else z<-z/totp # to ensure integrates to 1
  
  #  Get minimum number of squares to encompass pHR% probability
  nps<-length(p)
  zsort <- sort(as.vector(z)) 
  pt <- cumsum(zsort)/sum(zsort)
  
  
  #  p's for plotting and number of cells for calculating HR
  pplot<-matrix(0,nps,1)
  ncells<-matrix(0,nps,1)
  for(i in 1:nps){
    pplot[i]<-min(zsort[pt>=(1-p[i])])
    ncells[i]<-length(pt[pt>=(1-p[i])])
  }
  
  # Home range area
  HR<-dx*dy*ncells
  out<-list(HR, pplot, totp)
  names(out)<-c("HR", "ps", "totp")
  return(out)
}


### Function "Kernel_polyg_fast" to determine Kernel contours 
# Method: We propose here this function to identify Kernel contours. It re-uses part of R codes
#         proposed by Fieberg (2014) (http://dx.doi.org/10.13020/D6G59W), based on their work on
#         Utilisation Distribution (UD) by Fieberg and Kochanny (2005) (http://dx.doi.org/10.2193/0022-541X%282005%2969%5B1346:QHOTIO%5D2.0.CO;2).
#         The function Kernel_polyg_fast especially includes:
#             - the function "KernSmooth::bkde2D" to estimate kernel density
#             - the function calcHR to calculate kernel HR estimates for a particular contour
#             - different steps in spatial with packages "sp" and "sf" to create spatial polygons of kernel contours
#               (FYI: some steps are developped here to remove holes from polygons).
# Data input: a "sf" object of spatial locations [points] used to calculate Kernels.
# Parameters: - you can specify the grid_size (for "bkde2D" function) and the smoothing factor (H_value, in m if coordinates are projected).
#               (for details: https://www.rdocumentation.org/packages/KernSmooth/versions/2.23-18/topics/bkde2D )
#             - projection: you may want to add a projection for output, different from input (default: the same projection as in input is applied)
#             - prob: probabilities of Kernel contours (default: 95%, 75%, 50%).
# Packages needed: KernSmooth, sp, sf, grDevices, maptools.
# Function needed: calcHR.
# Output: spatial polygon of Kernel contours (in "sf" format).

Kernel_polyg_fast <- function(sf_points, grid_size=500, H_value=1500, prob=c(0.95, 0.75, 0.5), ID="no_id", projection=NULL){
  
  coords <- data.frame(sf::st_coordinates(sf_points)) # point coordinates
  
  if(is.null(projection)){ # define projection
    projection <- sf::st_crs(sf_points)
    projection <- projection$input
  }else{
    projection <- projection_selec
  }
  
  # Extrat min/max of lat/lon
  htemp<-apply(rbind(coords,coords),2,dpik)
  min.lon<-min(coords[,1])-20*htemp[1]
  max.lon<-max(coords[,1])+20*htemp[1]
  min.lat<-min(coords[,2])-20*htemp[2]
  max.lat<-max(coords[,2])+20*htemp[2]
  
  
  # Calculate Home Ranges with plug in method (separate bandwidths in x and y directions)
  UD1 <- KernSmooth::bkde2D(x=coords, bandwidth=H_value, gridsize=c(grid_size, grid_size), range.x=list(c(min.lon,max.lon),
                                                                                                        c(min.lat, max.lat)))
  
  # Determine kernel contours for each prob
  result_ALLprob <- NULL
  for (k in 1:length(prob)){
    
    prob_selec <- prob[k]
    
    # Extract XX% kernel contours using the function calcHR
    Coords_hr <- calcHR(UD1, p=prob_selec)
    
    # Get Percents_Kernel% conditional distributions
    dxdy <- (UD1$x1[2]-UD1$x1[1])*(UD1$x2[2]-UD1$x2[1])
    
    # Select UD
    UD1_selec <- UD1$fhat
    UD1_selec[UD1_selec<Coords_hr$ps[1]]<-0  
    UD1_selec <- UD1_selec/sum(UD1_selec*dxdy)
    
    # Extract contour lines of kernel
    contour_kernel1 <- grDevices::contourLines(UD1$x1, UD1$x2, UD1$fhat, levels=Coords_hr$ps[1,])
    
    # EXCLUDE HOLES in contours & do not keep contours with too much points
    # Different tests to clean the contour lines
    test <- NULL
    for(t in 1:length(contour_kernel1)){
      test <- c(test, length(contour_kernel1[[t]][[2]]))
    }
    
    if(sum(test >=4) == length(test)){ # If everything is ok
      
      shp_kernel1 <- maptools::ContourLines2SLDF(contour_kernel1)
      
      suppressWarnings(
        sp::proj4string(shp_kernel1) <- sp::CRS(projection)
      )
      
      # NEED sf PACKAGE: transform in sf object and calculate area
      SF1 <- sf::st_as_sf(shp_kernel1)
      SF1 <- sf::st_cast(SF1,"MULTILINESTRING")
      SF1 <- sf::st_cast(SF1,"MULTIPOLYGON")
      
      SF1 <- sf::st_make_valid(SF1)
      
      result_prob <- SF1
      
    }else{ # If it is not ok (e.g. a hole in the polygon)
      
      # Construct polygon manually (i = one line to define a polygon (or a hole !))
      for (i in 1:length(contour_kernel1)){
        
        contour_kernel1i <- list(contour_kernel1[[i]])
        shp_kernel1i <- maptools::ContourLines2SLDF(contour_kernel1i)
        
        if(length(contour_kernel1i[[1]]$x) >= 4){
          
          suppressWarnings(
            sp::proj4string(shp_kernel1i) <- sp::CRS(projection)
          )
          
          # NEED sf PACKAGE: transform in sf object and calculate area
          SFi <- sf::st_as_sf(shp_kernel1i)
          SFi <- sf::st_cast(SFi,"MULTILINESTRING")
          SFi <- sf::st_cast(SFi,"MULTIPOLYGON")
          
          
          if (i == 1){ # First polygon
            
            shp_kernel_final <- SFi
            
          }else{ # Next polygons
            
            if(isFALSE(exists("shp_kernel_final"))){
              shp_kernel_final <- SFi
            }else{
              
              Diff_test <- sf::st_intersection(shp_kernel_final, SFi)
              Union_test <- sf::st_union(shp_kernel_final, SFi)
              
              Area_test_shp <- ifelse(length(sf::st_area(shp_kernel_final))==0, 0, as.numeric(sf::st_area(shp_kernel_final)))
              Area_test_SFi <- ifelse(length(sf::st_area(SFi))==0, 0, as.numeric(sf::st_area(SFi)))
              Area_test_diff <- ifelse(length(sf::st_area(Diff_test))==0, 0, as.numeric(sf::st_area(Diff_test)))
              Area_test_union <- ifelse(length(sf::st_area(Union_test))==0, 0, sf::st_area(Union_test))
              
              # Remove hole OR union
              if (Area_test_diff == Area_test_SFi &
                  Area_test_SFi < Area_test_shp &
                  Area_test_union == Area_test_shp){ # If hole, if Area i+1 < Area i
                
                shp_kernel_final <- sf::st_sym_difference(shp_kernel_final, SFi)
                
              }else{
                
                if (Area_test_diff == Area_test_shp &
                    Area_test_SFi > Area_test_shp &
                    Area_test_union == Area_test_SFi){ # If hole, if Area i+1 > Area i
                  
                  shp_kernel_final <- sf::st_sym_difference(shp_kernel_final, SFi)
                  
                }else{
                  
                  shp_kernel_final <- sf::st_union(shp_kernel_final, SFi)
                  
                }
              }
            }
          }
        }else{}
      }
      
      shp_kernel_final <- shp_kernel_final[,1]
      result_prob <- shp_kernel_final
      
    }
    
    result_prob$level <- ID
    names(result_prob) <- c("ID", "geometry")
    result_prob$Prob <- prob_selec
    
    suppressWarnings(result_prob <-  # Resolve problem in sf format
                       sf::st_collection_extract(
                         result_prob, type = c("POLYGON", "POINT", "LINESTRING"),
                         warn = F))
    
    result_ALLprob <- rbind(result_ALLprob, result_prob)
    
  }
  
  return(result_ALLprob)
  
}


  ## 3 / ggplot settings to create a map with sf package

### Function "Set_limits_sf_ggplot" to determine X/Y min/max limits to create a map with a sf object in ggplot2
# You can set the expected limits:
#     - data_sf = a sf object
#     - X_prop & Y_prop = the proportions of map spaces expected in the plot. e.g. if you want a square, indicate X_prop = 1 and Y_prop = 1.
#     - percent_border = the percentage of "blank" space expected around the sf object. Default: 10%.
# Package needed: sf.
# Output: a simple data.frame of min/max of X and Y.

Set_limits_sf_ggplot <- function(data_sf, X_prop = 1.5, Y_prop = 1, percent_border = 10){
  
  X_rg <- c(min(sf::st_coordinates(st_union(data_sf))[,1]), max(sf::st_coordinates(st_union(data_sf))[,1]))
  Y_rg <- c(min(sf::st_coordinates(st_union(data_sf))[,2]), max(sf::st_coordinates(st_union(data_sf))[,2]))
  
  X_center <- mean(X_rg)
  Y_center <- mean(Y_rg)
  
  Range_X <- X_rg[2] - X_rg[1] 
  
  Range_Y <- Y_rg[2] - Y_rg[1]
  
  X_Y <- ((X_prop)/Y_prop)
  Y_X <- (Y_prop/(X_prop))
  
  if ((Range_X) / Range_Y >= X_Y){
    
    X_final <- c(X_center - ((Range_X/2)+((Range_X/2)*(percent_border/100))),
                 X_center + ((Range_X/2)+((Range_X/2)*(percent_border/100))))
    
    Y_final <- c(Y_center - (((Range_X/2)+((Range_X/2)*(percent_border/100)))*Y_X),
                 Y_center + (((Range_X/2)+((Range_X/2)*(percent_border/100)))*Y_X))
    
  }else{
    
    Y_final <- c(Y_center - ((Range_Y/2)+((Range_Y/2)*(percent_border/100))),
                 Y_center + ((Range_Y/2)+((Range_Y/2)*(percent_border/100))))
    
    X_final <- c(X_center - (((Range_Y/2)+((Range_Y/2)*(percent_border/100)))*X_Y),
                 X_center + (((Range_Y/2)+((Range_Y/2)*(percent_border/100)))*X_Y))
    
  }
  result <- matrix(c(X_final, Y_final), nrow = 2, ncol = 2)
  colnames(result) <- c("X", "Y")
  rownames(result) <- c("min", "max")
  result <- as.data.frame(result)
  return(result)
}

###########################################################################################################


### I // Data #############################################################################################
  ## 1 / Direction 
  Direction <- ".../Planque_et_al_Foraging_areas_V1_2"
  #Direction <- "C:/Users/yplanq01/Documents/CEBC/Article_Planque_et_al_Niche_Overlap/Scripts/Foraging_areas/Planque_et_al_Foraging_areas_V1_2"

  ## 2 / Import data
  #  Data used in this script are available on SEANOE: 
  #   Planque Yann, Caurant Florence, Vincent Cecile (2021).
  #   Dive data obtained from telemetry tracking of ten harbour seals (Phoca vitulina)
  #   and twelve grey seals (Halichoerus grypus), captured in the Baie de Somme, France,
  #   in 2008 and 2012, and fitted with GPS/GSM tags. SEANOE.
  #   https://doi.org/10.17882/80016
  
  # >> Please place the data (cf. 'Dive data' 83096.csv in SEANOE) in the subfolder "Input" of your direction

  #   Pv : Phoca vitulina (harbour seals)  //  Hg : Halichoerus grypus (grey seals)
  Dive_data <- read.table(paste0(Direction, "/Input/","83096.csv"),dec=".", sep=";", header = T)
  head(Dive_data)
  
  # Here, remove two individuals before the analyses:
      # G09 (tag malfunction after 5 days) and S09 (low tracking duration, and low resolution at its scale)
      Dive_data <- Dive_data[!Dive_data$Seal_ID %in% c("G09", "S09"),]
  
  Dive_data$ds_date <- as.POSIXct(Dive_data$ds_date, format="%Y-%m-%d %H:%M:%S",tz="UTC")
  
  # Seal individuals
  Seal_IDs <- unique(Dive_data$Seal_ID)
  Seal_IDs
  
###########################################################################################################

  
### II // Identification of likely foraging dives of harbour and grey seals ###############################
## Method: cf. 'vertical approach' presented by Planque et al. (2020) (https://doi.org/10.1007/s00227-019-3636-8)
## Objective: selection of the faster U-shaped dives (~22.5% by individual)
## Note: the method used here aimed at selecting faster U-shaped dives, that are expected to
##       describe the benthic foraging behaviour of harbour and grey seals (cf. Planque et al., 2020). 
##       In the presented case study, we are confident in stating that both species should exhibit
##       this behaviour as their diet are mostly composed of benthic flatfish
##       (Spitz et al., 2015 https://doi.org/10.1051/alr/2015001; Planque et al., in review).
##       Thus, such an approach may only be suitable for benthic diving predators. 
  
  ## 1 / Dive shape 
    # A / Calculate MCTS and TAD at the individual level
  
  Dive_data2 <- NULL
    for (i in 1:length(Seal_IDs)){ # For each individual...
      
      Dive_data_ind <- Dive_data[Dive_data$Seal_ID == Seal_IDs[i],]
      
      # ... calculate the MCTS
      Dive_data_ind$MCTS <- 
      MCTS_ind(Dive_data_ind, DIVE_DUR = "DIVE_DUR", MAX_DEP = "MAX_DEP", PERCENT_AREA = "PERCENT_AREA", 
                MCTS_min=1.0, MCTS_max=3.0, MCTS_by=0.1, Dive_threshold=1.5,
                return_plot_obj=F, show_plot=F)
      
      # Create MCTS plot and export
      # Plot_MCTS <-
      # MCTS_ind(Dive_data_ind, DIVE_DUR = "DIVE_DUR", MAX_DEP = "MAX_DEP", PERCENT_AREA = "PERCENT_AREA",
      #          MCTS_min=1.0, MCTS_max=3.0, MCTS_by=0.1, Dive_threshold=1.5,
      #          return_plot_obj=T, show_plot=F)
      # 
      # ggsave(Plot_MCTS,
      #        filename = paste0(Direction, "/Plot/MCTS/", "00_MCTS_", Seal_IDs[i],".png"), dpi = 300,
      #        width = 11, height = 7)
      
      
      # ... calculate the TAD for each dive
      Dive_data_ind$TAD <- 
        TAD_index(Dive_data_ind, DIVE_DUR = "DIVE_DUR", MAX_DEP = "MAX_DEP",
                 PERCENT_AREA = "PERCENT_AREA", MCTS = "MCTS", Dive_threshold=1.5)
      
      Dive_data2 <- rbind(Dive_data2, Dive_data_ind)  # New dive data with individual MCTS and TAD
      
    }; rm(Dive_data_ind, Plot_MCTS)
    
  
    # B / Identify the most U-shaped dives 
    # We define here TAD thresholds to select ~25% of the most U-shaped dives for each individual

    # Define quantile at 75% to determine TAD thresholds
    TAD_thres <- as.data.frame(
                  Dive_data2 %>% dplyr::group_by(Seal_ID) %>%
                    dplyr::summarise(TAD_thres = quantile(TAD, 0.75)) %>%
                    dplyr::mutate(TAD_thres = as.numeric(TAD_thres)))
    
    
    # Add this to the data.frame
    Dive_data2 <- dplyr::left_join(Dive_data2, TAD_thres)
    
    # U-dives column
    Dive_data2$U_dives <- ifelse(Dive_data2$TAD >= Dive_data2$TAD_thres, 1, 0)
  
    # Check if it is 25% of all dives
    sum(Dive_data2$U_dives) / nrow(Dive_data2)
    
    
    
  ## 2 / Vertical descent speed 
    # A / V_DESC : vertical descent speed during the first T% of dive duration 
    #     Calculated between the surface (here 1.5 m) and intermediate depth D1 (at 10% of dive duration)
    #     in m.s-1
    
    Dive_data2$V_DESC <- (Dive_data2$D1 - 1.5) / (Dive_data2$DIVE_DUR * (10/100))
  
    # can be negative when D1 is above 1.5 m (for very short and shallow dives).
    # In this case, V_DESC is 0. Clean that:
    Dive_data2$V_DESC <- ifelse(Dive_data2$V_DESC < 0, 0, Dive_data2$V_DESC)
    
    
    # B / Identify the faster U-shaped dives
    # Objective: exclude the slowest U-shaped that are likely associated with a resting/sleeping behaviour
    # Method: we define here V_DESC thresholds to exclude 10% of the previously selected U-shaped dives,
    #         associated to the lowest vertical descent speeds
    
    # Define quantile at 10% to determine V_DESC thresholds (only on U-shaped dives)
    V_DESC_thres <- as.data.frame(
      Dive_data2 %>% dplyr::filter(U_dives == 1) %>% dplyr::group_by(Seal_ID) %>%
        dplyr::summarise(V_DESC_thres = quantile(V_DESC, 0.1)) %>%
        dplyr::mutate(V_DESC_thres = as.numeric(V_DESC_thres)))
    
    # Add this to the data.frame
    Dive_data2 <- dplyr::left_join(Dive_data2, V_DESC_thres)
    
    # Fast-U-dives column
    Dive_data2$Fast_U_dives <- ifelse(Dive_data2$TAD >= Dive_data2$TAD_thres &
                                      Dive_data2$V_DESC >= Dive_data2$V_DESC_thres, 1, 0)
    
    # Check if it is ~ 22.5% of all dives
    sum(Dive_data2$Fast_U_dives) / nrow(Dive_data2)

# Export a table of thresholds  
Out_summary <- as.data.frame(
  Dive_data2 %>% dplyr::group_by(Species, Seal_ID) %>%
    dplyr::summarise(N_dives = n(), MTCS = mean(MCTS),
                     TAD_thres = mean(TAD_thres), V_DESC_thres = mean(V_DESC_thres),
                     N_fast_U_dives = sum(Fast_U_dives)))

write.table(Out_summary, paste0(Direction, "/Output/Dive_criteria_summary.csv"), dec=".", sep=";", row.names = F)
    
###########################################################################################################
    
  
### III // Identification of likely foraging areas ########################################################
## => Spatialise the selected likely foraging dives (cf. II) with Kernel contours
## Method: See details on the function "Kernel_polyg_fast" in section the "0 // Useful functions" above.
##         It is partly based on R codes provided by Fieberg (2014) (http://dx.doi.org/10.13020/D6G59W).
##         New codes are implemented to extract spatial polygons in "sf" format (usefull for ggplot).
## Objective: calculate Kernel contours at 95%, 75% and 50% of likely foraging dives, and create maps at species/individual levels.


  ## 1 / Prepare dive data
  # Select likely foraging dives, and only keep necessary columns
  Dive_data_foraging <- Dive_data2 %>% dplyr::filter(Fast_U_dives == 1) %>%
          dplyr::select(FID_dive, Species, Seal_ID, Trip_nb, ds_date, start_lat, start_lon) 


  ## 2 / Define a spatial projection
  # Many projections are available on this website: http://spatialreference.org/ref/
  #     => click on "Proj4" to obtain the code of the chosen projection

  # Study projection, here:  RGF 1993 Lambert 93 (ESRI format) [http://spatialreference.org/ref/esri/102110/]
  projection_selec <- "+proj=lcc +lat_1=44 +lat_2=49 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +units=m +no_defs" 


  ## 3 / Base map for this study case 
  # You can find base map here: https://ec.europa.eu/eurostat/fr/web/gisco/geodata/reference-data/administrative-units-statistical-units/countries
  # Here we dowloaded the shapefile CNTR_BN_01M_2020_3035.shp (countries 2020 => scale 1:1 Million => SHP)
  # Place your shapefile in the subfolder ".../Input/Base_map"
  # Now read your shapefile...
  Base_map <- rgdal::readOGR(dsn=paste(Direction,"/Input/Base_map",sep=""), layer="CNTR_RG_01M_2020_4326")
  
  # Select countries in our study area:
  EU_cntry <- c("AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "IRL",
                "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE",
                "GBR", "CHE", "NOR") # and welcoming non-EU countries!
  Europe <- Base_map[Base_map$ISO3_CODE %in% EU_cntry,]
  
        # Another basemap...
        # Europe_lbrt93 <- rgdal::readOGR(dsn=paste(Direction,"/Input/Base_map",sep=""),layer="Fond_Europe_LAMBERT1993_selection_pays")

  Europe <- sp::spTransform(Europe, sp::CRS(projection_selec)) # project the base map
  Europe_proj_sf <- sf::st_as_sf(Europe) # convert it in "sf" object
  
  ggplot() + # Check the base map
    geom_sf(data=Europe_proj_sf, fill="#D9D9D9", color="#0B0B0B", size=0.4) +
    theme_bw(base_size = 14)

  
  ## 4 / Project foraging dives 
  sp::coordinates(Dive_data_foraging) <- c("start_lon", "start_lat")
  sp::proj4string(Dive_data_foraging) <- sp::CRS("+proj=longlat +datum=WGS84") # project in WGS84 (cf. coordinates in decimal degrees)
  Dive_data_foraging <- sp::spTransform(Dive_data_foraging, CRS(projection_selec)) # transform with our projection
  Dive_data_foraging_sf <- sf::st_as_sf(Dive_data_foraging)
  
    
  # Have a look to the dive points... 

        # Select part of dives [EXAMPLE for grey seals Hg]
        Dive_to_plot <- Dive_data_foraging_sf[Dive_data_foraging_sf$Species == "Hg",]
        Dive_to_plot <- Dive_to_plot %>% arrange(Seal_ID, ds_date) 
        
        # Unique trip_nb bys inds.
        Dive_to_plot$Trip_nb_Ind <- paste0(Dive_to_plot$Seal_ID, "_", Dive_to_plot$Trip_nb)
        
        # Map limits & scalevar size
        Lims_dives <- Set_limits_sf_ggplot(data_sf = Dive_to_plot, X_prop = 1.5, Y_prop = 1, percent_border = 10)
        Size_scalebar_dives <- round(((abs(Lims_dives$X[2] - Lims_dives$X[1])/1000)/6)/20)*20  # round at the closest 20 km of 1/6 of X range  # auto x2 in scalebar()
        
        # The map
        Plot_dives <- 
        ggplot() +
          # Dive points
          geom_sf(data=Dive_to_plot,
                  aes(color=Seal_ID), shape=16, size=0.55, alpha=0.1) +
          # Base map
          geom_sf(data = Europe_proj_sf, fill="#D9D9D9", color="#0B0B0B", size=0.4) +
          # We could also looked at trips with geom_path. But it is more accurate to do so with all locations (i.e. without dive selection) 
          # geom_path(data = cbind(st_coordinates(Dive_to_plot),
          #                        data.frame(Seal_ID=Dive_to_plot$Seal_ID,
          #                                   Trip_nb_Ind=Dive_to_plot$Trip_nb_Ind)),
          #           aes(x=X, y=Y, color=Seal_ID, group=Trip_nb_Ind), size=0.1) +
          
          # Set colors
          scale_color_manual(name="Inds.",                                                   
                             values=rainbow(length(unique(Dive_to_plot$Seal_ID)))) +
          # Map limits
          coord_sf(crs = projection_selec,
                   xlim = Lims_dives$X,
                   ylim = Lims_dives$Y,
                   expand = FALSE) +
          # Scale bar
          ggsn::scalebar(x.min = Lims_dives$X[1],
                         x.max = Lims_dives$X[2] - (Lims_dives$X[2] - Lims_dives$X[1])*0.05,
                         y.min = Lims_dives$Y[1] + (Lims_dives$Y[2] - Lims_dives$Y[1])*0.0625,
                         y.max = Lims_dives$Y[2],
                         dist = Size_scalebar_dives, dist_unit = "km", st.dist = 0.0325,
                         st.size=3.25, border.size = 0.25,
                         transform = FALSE) +
          # Other plot parameters/themes
          labs(x = "", y = "", title = paste0("Likely foraging dives of grey seals (N = ", nrow(Dive_to_plot),")")) +
          theme_bw(base_size = 14) +
          theme(axis.text.x = element_text(size = 11, color="black"),
                axis.text.y = element_text(size = 11, color="black"),
                axis.title = element_text(size = 13, color="black"),
                plot.title = element_text(size = 14, color="black", face = "bold"),
                legend.title = element_text(size = 13, color="black", face = "bold"),
                legend.text = element_text(size = 13, color="black")) +
          guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) + # Set points in legends
          geom_blank()
        
        # Export Map
        ggsave(Plot_dives,
               filename = paste0(Direction, "/Plot/", "00_Dives_Hg.png"), dpi = 300,
               width = 7.5, height = 4.5)
  
  
  ## 5 / Foraging areas (Kernels at 95%, 75%, 50%)
  # Settings for Kernels
  # Different smoothing factors (H) between harbour (Pv) and grey seals (Hg), to consider their differences in spatial extent.
  
  grid_size <- 500
  H_value_Pv <- 1000
  H_value_Hg <- 2500
  
  Probs_test <- c(0.95, 0.75, 0.5)
  
  COLOR_Pv <- "#4DE600"
  COLOR_Hg <- "#0000DE"
  
  
    # A / For all individuals pooled (by species)
  
    # Harbour seals_________
    # Kernels
    Kernels_all_Pv <- 
      Kernel_polyg_fast(sf_points = Dive_data_foraging_sf[Dive_data_foraging_sf$Species == "Pv",],
                        ID = "Pv", grid_size = grid_size, H_value = H_value_Pv, prob = Probs_test)
  
    # Map limits & scalevar size
    Lims_all_Pv <- Set_limits_sf_ggplot(data_sf = Kernels_all_Pv, X_prop = 1.5, Y_prop = 1, percent_border = 10)
    Size_scalebar_Pv <- round(((abs(Lims_all_Pv$X[2] - Lims_all_Pv$X[1])/1000)/8)/20)*20  # round at the closest 20 km of 1/8 of X range  # auto x2 in scalebar()
    
    # Plot that...
    Map_all_Pv <- 
    ggplot() +
      # Kernels
      geom_sf(data=Kernels_all_Pv,
              aes(alpha=factor(Prob, levels = Probs_test)),
              fill=COLOR_Pv, color=COLOR_Pv, size=0.4) +
      # Base map
      geom_sf(data = Europe_proj_sf, fill="#D9D9D9", color="#0B0B0B", size=0.4) +
      
      # Transparency of kernels
      scale_alpha_manual(name="Kernel\ndensities",                                                   
                         values=c(0.2, 0.5, 1), labels = paste0(Probs_test*100,"%")) +
      # Map limits
      coord_sf(crs = projection_selec,
               xlim = Lims_all_Pv$X,
               ylim = Lims_all_Pv$Y,
               expand = FALSE) +
      # Scale bar
      ggsn::scalebar(x.min = Lims_all_Pv$X[1],
                     x.max = Lims_all_Pv$X[2] - (Lims_all_Pv$X[2] - Lims_all_Pv$X[1])*0.05,
                     y.min = Lims_all_Pv$Y[1] + (Lims_all_Pv$Y[2] - Lims_all_Pv$Y[1])*0.0625,
                     y.max = Lims_all_Pv$Y[2],
                     dist = Size_scalebar_Pv, dist_unit = "km", st.dist = 0.0325,
                     st.size=3.25, border.size = 0.25,
                     transform = FALSE) +
      # Other plot parameters/themes
      labs(x = "", y = "", title = paste0(" All harbour seals (", length(unique(
             Dive_data_foraging_sf[Dive_data_foraging_sf$Species == "Pv",]$Seal_ID)),
             " individuals)")) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(size = 11, color="black"),
            axis.text.y = element_text(size = 11, color="black"),
            axis.title = element_text(size = 13, color="black"),
            plot.title = element_text(size = 14, color="black", face = "bold"),
            legend.title = element_text(size = 13, color="black", face = "bold"),
            legend.text = element_text(size = 13, color="black")) +
      geom_blank()
    
    Map_all_Pv
    
    # Export Map
    ggsave(Map_all_Pv,
           filename = paste0(Direction, "/Plot/", "01_Kernels_all_Pv.png"), dpi = 300,
           width = 7.5, height = 4.5)
    
    # Save polygons in a shapefile
    sf::st_write(Kernels_all_Pv, paste0(Direction, "/Output/Shapefiles/01_Kernels_all_Pv.shp"))
    
    
    
    
    # Grey seals_________
    # Kernels
    Kernels_all_Hg <- 
      Kernel_polyg_fast(sf_points = Dive_data_foraging_sf[Dive_data_foraging_sf$Species == "Hg",],
                        ID = "Hg", grid_size = grid_size, H_value = H_value_Hg, prob = Probs_test)
    
    # Map limits & scalevar size
    Lims_all_Hg <- Set_limits_sf_ggplot(data_sf = Kernels_all_Hg, X_prop = 1.5, Y_prop = 1, percent_border = 10)
    Size_scalebar_Hg <- round(((abs(Lims_all_Hg$X[2] - Lims_all_Hg$X[1])/1000)/8)/20)*20  # round at the closest 20 km of 1/8 of X range  # auto x2 in scalebar()
    
    # Plot that...
    Map_all_Hg <- 
      ggplot() +
      geom_sf(data=Kernels_all_Hg,
              aes(alpha=factor(Prob, levels = Probs_test)),
              fill=COLOR_Hg, color=COLOR_Hg, size=0.4) +
      geom_sf(data = Europe_proj_sf, fill="#D9D9D9", color="#0B0B0B", size=0.4) +
      scale_alpha_manual(name="Kernel\ndensities",                                                   
                         values=c(0.2, 0.5, 1), labels = paste0(Probs_test*100,"%")) +
      coord_sf(crs = projection_selec,
               xlim = Lims_all_Hg$X,
               ylim = Lims_all_Hg$Y,
               expand = FALSE) +
      ggsn::scalebar(x.min = Lims_all_Hg$X[1],
                     x.max = Lims_all_Hg$X[2] - (Lims_all_Hg$X[2] - Lims_all_Hg$X[1])*0.05,
                     y.min = Lims_all_Hg$Y[1] + (Lims_all_Hg$Y[2] - Lims_all_Hg$Y[1])*0.0625,
                     y.max = Lims_all_Hg$Y[2],
                     dist = Size_scalebar_Hg, dist_unit = "km", st.dist = 0.0325,
                     st.size=3.25, border.size = 0.25,
                     transform = FALSE) +
      labs(x = "", y = "", title = paste0(" Grey seals (", length(unique(
        Dive_data_foraging_sf[Dive_data_foraging_sf$Species == "Hg",]$Seal_ID)),
        " individuals)")) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(size = 11, color="black"),
            axis.text.y = element_text(size = 11, color="black"),
            axis.title = element_text(size = 13, color="black"),
            plot.title = element_text(size = 14, color="black", face = "bold"),
            legend.title = element_text(size = 13, color="black", face = "bold"),
            legend.text = element_text(size = 13, color="black")) +
      geom_blank()
    
    Map_all_Hg
    
    # Export Map
    ggsave(Map_all_Hg,
           filename = paste0(Direction, "/Plot/", "01_Kernels_all_Hg.png"), dpi = 300,
           width = 7.5, height = 4.5)
    
    # Save polygons in a shapefile
    sf::st_write(Kernels_all_Hg, paste0(Direction, "/Output/Shapefiles/01_Kernels_all_Hg.shp"))

    
    # B / At the individual level
  
    Seal_IDs # all inds
    
    Kernels_all_inds <- NULL
    for (u in 1:length(Seal_IDs)){  # individual loop
      
      # Individual data & settings
      Dive_data_sf_ind <- Dive_data_foraging_sf[Dive_data_foraging_sf$Seal_ID == Seal_IDs[u],]
      
      if(unique(Dive_data_sf_ind$Species) == "Pv"){
        SPP <- "Harbour seal"
        COLOR_ind <- COLOR_Pv
        H_value_ind <- H_value_Pv}
      
      if(unique(Dive_data_sf_ind$Species) == "Hg"){
        SPP <- "Grey seal"
        COLOR_ind <- COLOR_Hg
        H_value_ind <- H_value_Hg}
      
      
      # Kernels
      Kernels_ind <- 
        Kernel_polyg_fast(sf_points = Dive_data_sf_ind,
                          ID = Seal_IDs[u], grid_size = grid_size, H_value = H_value_ind, prob = Probs_test)
      
      Kernels_ind$Species <- unique(Dive_data_sf_ind$Species)
      
      
      # Map limits & scalevar size
      Lims_ind <- Set_limits_sf_ggplot(data_sf = Kernels_ind, X_prop = 1.5, Y_prop = 1, percent_border = 10)
      
      if(unique(Dive_data_sf_ind$Species) == "Pv"){
        Size_scalebar_ind <- ceiling(((abs(Lims_ind$X[2] - Lims_ind$X[1])/1000)/8)/5)*5 }  # round at the closest 5 km of 1/8 of X range  # auto x2 in scalebar()
      
      if(unique(Dive_data_sf_ind$Species) == "Hg"){
        Size_scalebar_ind <- ceiling(((abs(Lims_ind$X[2] - Lims_ind$X[1])/1000)/8)/10)*10 }  # round at the closest 10 km of 1/8 of X range  # auto x2 in scalebar()
      
       
      # Plot that...
      Map_ind <- 
        ggplot() +
        geom_sf(data=Kernels_ind,
                aes(alpha=factor(Prob, levels = Probs_test)),
                fill=COLOR_ind, color=COLOR_ind, size=0.4) +
        geom_sf(data = Europe_proj_sf, fill="#D9D9D9", color="#0B0B0B", size=0.4) +
        scale_alpha_manual(name="Kernel\ndensities",                                                   
                           values=c(0.2, 0.5, 1), labels = paste0(Probs_test*100,"%")) +
        coord_sf(crs = projection_selec,
                 xlim = Lims_ind$X,
                 ylim = Lims_ind$Y,
                 expand = FALSE) +
        ggsn::scalebar(x.min = Lims_ind$X[1],
                       x.max = Lims_ind$X[2] - (Lims_ind$X[2] - Lims_ind$X[1])*0.05,
                       y.min = Lims_ind$Y[1] + (Lims_ind$Y[2] - Lims_ind$Y[1])*0.0625,
                       y.max = Lims_ind$Y[2],
                       dist = Size_scalebar_ind, dist_unit = "km", st.dist = 0.0325,
                       st.size=3.25, border.size = 0.25,
                       transform = FALSE) +
        labs(x = "", y = "", title = paste0(SPP, " ", Seal_IDs[u])) +
        theme_bw(base_size = 14) +
        theme(axis.text.x = element_text(size = 11, color="black"),
              axis.text.y = element_text(size = 11, color="black"),
              axis.title = element_text(size = 13, color="black"),
              plot.title = element_text(size = 14, color="black", face = "bold"),
              legend.title = element_text(size = 13, color="black", face = "bold"),
              legend.text = element_text(size = 13, color="black")) +
        geom_blank()
      
      Map_ind
      
      # Export Map
      ggsave(Map_ind,
             filename = paste0(Direction, "/Plot/", "02_Kernels_", Seal_IDs[u],".png"), dpi = 300,
             width = 7.5, height = 4.5)
      
      # Merge ind. data together    
      Kernels_all_inds <- rbind(Kernels_all_inds, Kernels_ind)
      
      rm(Kernels_ind, Dive_data_sf_ind, SPP, COLOR_ind, H_value_ind)
      
    }; rm(u, i) # End of individual loop
    
    
    # Save polygons in a shapefile
    sf::st_write(Kernels_all_inds[Kernels_all_inds$Species == "Pv",], paste0(Direction, "/Output/Shapefiles/02_Kernels_Pv_all_inds.shp"))
    sf::st_write(Kernels_all_inds[Kernels_all_inds$Species == "Hg",], paste0(Direction, "/Output/Shapefiles/02_Kernels_Hg_all_inds.shp"))
    
    

# Save the data in our session
save.image(paste(Direction, "Output", "Foraging_areas.RData", sep = "/"), safe = TRUE,
           compress = "xz")
#load(paste(Direction, "Output", "Foraging_areas.RData", sep = "/"))


# Export the analysed dive data in a table
write.table(Dive_data2, paste0(Direction, "/Output/Seal_Dives_BDS_data_analysed.csv"), dec=".", sep=";", row.names = F)

########################################################################################################### 


### IV // OPTIONAL: Visualise intermediate depths (D1:D9) with a time-depth profile #######################
#  Show time-depth profile for one individual during a specified period
#  e.g. grey seal G08, 2012-09-22 from 00:00 to
Time_depth_example <- 
  Dive_data2 %>% filter(Seal_ID == "G08" &
                         # Choose dates of a period at-sea
                         ds_date >= as.POSIXct("2012-09-22 00:00:00", format="%Y-%m-%d %H:%M:%S",tz="UTC") &
                         ds_date <= as.POSIXct("2012-09-22 01:00:00", format="%Y-%m-%d %H:%M:%S",tz="UTC")) %>%
  mutate(D0=1.5, D10=1.5) %>%  # add surface points
  # gather to have intermediate depths in long data.frame
  gather(key = "D_depth", value = "Depth_m", D1:D9, D0, D10) %>%  # select intermediate depths
  mutate(Percent_Depth = as.numeric(substring(D_depth, 2))*10) %>%  # % of Dive duration, here D1:D9 by 10%. D0 at 0% and D10 at 100%.
  mutate(D_date = ds_date + (DIVE_DUR*Percent_Depth)/100) %>%
  arrange(Seal_ID, FID_dive, D_date) %>%
  ggplot(aes(x=D_date, y=-Depth_m)) +
  geom_hline(yintercept=0, color = "#0053CC", size=0.35) +  # surface
  geom_hline(yintercept=-1.5, color = "#0053CC", linetype="dashed", size=0.5) +  # dive detection (1.5 m)
  geom_line() +
  geom_point(fill="red2", shape=21, size=2) +
  labs(x = "Time", y = "Depth (m)") +  # title = "Grey seal G08", subtitle = "Time-depth profile 2012-09-22 from 00:00 to 01:00"
  scale_y_continuous(labels = abs) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(size = 11, color="black"),
        axis.text.y = element_text(size = 11, color="black"),
        axis.title = element_text(size = 13, color="black"),
        plot.title = element_text(size = 14, color="black", face = "bold"),
        legend.title = element_text(size = 13, color="black", face = "bold"),
        legend.text = element_text(size = 13, color="black"),
        panel.grid = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5)))

Time_depth_example

ggsave(Time_depth_example,
       filename = paste0(Direction, "/Plot/", "03_SUPP_Time_depth_ex.png"), dpi = 300,
       width = 9, height = 3.5)

###########################################################################################################


### Additional information ################################################################################
##
## >> The individual foraging areas presented in Planque et al. (2021) were identified using this
##    script. Only individuals for which we had a whisker and measured stable isotopes along it, were
##    presented in this study. i.e. 8 harbour seals and 10 grey seals. Note that foraging areas of the
##    grey seal G09 could not be identified due to a tag malfunction (only 5 days after the start tracking).
##
## >> The maps in Planque et al. (2021) were realised with ArcMap 10.6, using the exported shapefiles.
##
## >> The foraging areas presented in Planque et al. (2020) were identified in a specific study area
##    (Eastern English Channel), during a specific period (excluding breeding & moulting periods),
##    due to methodological choices for comparing two approaches. Thus, results in this script may vary
##    to those in this article. This script therefore re-used the vertical approach presented in this article.
##
########################################################################################################### 


###########################################################################################################
##################################### END #################################################################
###########################################################################################################