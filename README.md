# Identification of foraging areas with dive data obtained from telemetry tracking
## Selection of likely foraging dives with a vertical approach [(Planque et al., 2020)](https://doi.org/10.1007/s00227-019-3636-8)
• **OBJECTIVES** : Identify the foraging areas of diving predators (here harbour and grey seals) by analysing dive data obtained from telemetry tags (here GPS/GSM tags). The two species are usually considered as benthic foragers and therefore typically perform fast U-shaped dives to forage. We therefore aimed at selecting this type of dives (called 'likely foraging dives') using two dive criteria

• **METHODS** : 
<br>(1) Calculate individual **Minimum Cost of Transport Speed (MCTS)**. It is based on TAD index calculation, and should ideally be run at the individual level (method presented by [Vincent et al., 2016](https://doi.org/10.1093/icesjms/fsw102)). Function *MCTS_ind*. </br>
<br>(2) Determine dive shape with the **Time Allocation at Depth (TAD)** index [(Fedak et al., 2001)](https://doi.org/10.1111/j.1748-7692.2001.tb00982.x). Function *TAD_index*.</br>
<br>(3) Select **likely foraging dives** using two dive criteria: dive shape (TAD) and vertical descent speed (vertical approach by [Planque et al., 2020](https://doi.org/10.1007/s00227-019-3636-8)).</br>
<br>(4) **Spatialise foraging areas** with likely foraging dives, by determining **Kernel density contours** at 95%, 75% and 50%. Kernels were identified here with the function *Kernel_polyg_fast*, based on some part of codes proposed by [Fieberg (2014)](http://dx.doi.org/10.13020/D6G59W), and on new code using functions in "sf" package. Output: "sf" polygons.</br>

• **LANGUAGE** : R (version 4.0.2).

• **CASE STUDY** : 9 harbour seals and 11 grey seals captured in the *Baie de Somme* (France), respectively in 2008 and 2012, and fitted with GPS/GSM tags.

## Script developped as part of: 
Planque, Y., Spitz, J., Authier, M., Vincent, C., Caurant, F. Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus) at their Southern European limit range (Eastern English Channel). Submitted in 'Ecology and Evolution'.

Preprint version (1st version submitted in 2020-11): https://doi.org/10.22541/au.160508195.50224560/v1

Peer review status (2021-03-31): minor revision

## Data used in this script are freely available on SEANOE:
Planque Yann, Caurant Florence, Vincent Cécile (2021). Dive data obtained from telemetry tracking of ten harbour seals (Phoca vitulina) and twelve grey seals (Halichoerus grypus), captured in the Baie de Somme, France, in 2008 and 2012, and fitted with GPS/GSM tags. SEANOE. https://doi.org/10.17882/80016

Before running the script, please download the file repository in the ZIP file and place it on your desktop. Place the data previously dowloaded in the subfolder "Input".

First publication on GitHub : 2020-11-05

Last update : 2021-04-12 (Version 1.3)

## Author : Yann Planque(1)
 Affiliation :
    (1) Centre d'Etudes Biologiques de Chizé (CEBC, UMR 7372 CNRS - La Rochelle Université), La Rochelle, France

## Contact : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
