# ravis_1006_final_project

This repo contains the final project data, code, and deliverables for Gov 1006. 

The paper to be replicated is "Titling community land to prevent deforestation: An evaluation of a best-case program in Morona-Santiago, Ecuador" by Butaine, et al, 2015. A PDF of the paper can be found in the top level directory in the repo. The original code and data is all located in the "TitlingDeforestation" folder.

This is an evaluation of the impact of a donor-funded land titling and land management program for indigenous communities in Ecuador. They believe that this is virtually a best-case policy intervention. The data for the model was generated from spatial data on intervention areas and title boundaries provided by USAID and NGOs and publicly available spatial raster data such as population density (Landscan) and forest cover (Global Forest Change (GFC) by Hansen et al). The change in GFC in an area is used as the metric for deforestation. Areas receiving the program treatment were matched with areas that did not receive the treatment using a genetic matching algorithm. The treatment effect over the five years after is then estimated with a difference-in-difference OLS model. They find that treatment does not have a significant effect on deforestation rates

The RMDs for the milestones are named "Milestone X."

The commented code for Milestone 3 is in "Rep Code Analysis.R"
