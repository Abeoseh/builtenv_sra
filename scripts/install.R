if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib="../my_R_libs", repos="https://cloud.r-project.org")

BiocManager::install("biomformat",lib="../my_R_libs")
