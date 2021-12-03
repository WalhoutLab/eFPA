#install.packages("devtools")
require(devtools)

metanr_packages()

# Step 2: Install MetaboAnalystR without documentation
devtools::install_github("xia-lab/MetaboAnalystR", build = FALSE, build_vignettes = FALSE)






install.packages('Cairo')
library(Cairo)

library(MetaboAnalystR)
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine")
mSet<-InitDataObjects("conc", "msetora", FALSE)
mSet<-Setup.MapData(mSet, tmp.vec);
mSet<-CrossReferencing(mSet, "name");

mSet<-SetCurrentMsetLib(mSet, "location", 0)

#==> same as what we used before 