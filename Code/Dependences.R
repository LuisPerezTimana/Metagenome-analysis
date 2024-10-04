##### MICROECO

## DEPENDENCES

# CRAN

# allow more waiting time to download each package
options(timeout = 1000)
# If a package is not installed, it will be installed from CRAN
# First select the packages of interest
tmp <- c("microeco", "mecoturn", "MASS", "GUniFrac", "ggpubr", "randomForest", "ggdendro", "ggrepel", "agricolae", "igraph", "picante", "pheatmap", "rgexf", 
         "ggalluvial", "ggh4x", "rcompanion", "FSA", "gridExtra", "aplot", "NST", "GGally", "ggraph", "networkD3", "poweRlaw", "ggtern", "SRS", "performance")
# Now check or install
for(x in tmp){
  if(!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
}

# BIOCONDUCTOR

install.packages("BiocManager")
install.packages("file2meco", repos = BiocManager::repositories())
install.packages("MicrobiomeStat", repos = BiocManager::repositories())
install.packages("WGCNA", repos = BiocManager::repositories())
BiocManager::install("ggtree")
BiocManager::install("metagenomeSeq")
BiocManager::install("ALDEx2")
BiocManager::install("ANCOMBC")

# GITHUB
# download link of the compressed packages archive
# Alternative from Gitee "https://gitee.com/chiliubio/microeco_dependence/releases/download/v0.20.0/microeco_dependence.zip"
url <- "https://github.com/ChiLiubio/microeco_dependence/releases/download/v0.20.0/microeco_dependence.zip"
# allow more time to download the zip file in R
options(timeout = 2000)
# Another way is to open the upper url in browser to download the zip file and move it to the current R working directory
download.file(url = url, destfile = "microeco_dependence.zip")
# uncompress the file in R
tmp <- "microeco_dependence"
unzip(paste0(tmp, ".zip"))
# install devtools
if(!require("devtools", character.only = TRUE)){install.packages("devtools", dependencies = TRUE)}
# run these one by one
devtools::install_local(paste0(tmp, "/", "SpiecEasi-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "mixedCCA-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "SPRING-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "NetCoMi-main.zip"), repos = BiocManager::repositories())
devtools::install_local(paste0(tmp, "/", "beem-static-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "chorddiag-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggradar-master.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggnested-main.zip"), dependencies = TRUE)
devtools::install_local(paste0(tmp, "/", "ggcor-1-master.zip"), dependencies = TRUE)

