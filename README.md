#TCGAGeneReport

This is for educational and research purposes only. 

http://bioinformatics.mdanderson.org

The goal of the Gene Survey project is to take TCGA Standardized Data and make gene level cross-disease and sometimes cross platform comparisons and evaluations simpler to perform, as these represent common questions asked of PIs and analysts. That is, we are asnwering a common question: how does gene XYZ act across disease types?

This is the R package used for this project. The Java package is available at https://github.com/GeneSurvey/TcgaGSData but the JAR is included with this source.

The data file will be available soon.

#General Pre-Requisites

The TCGAGeneReport package requires Java 1.8.

These instructions have been testing using R 3.2.3 and R 3.2.2.

#Linux Prerequisites

Linux (Debian/Ubuntu) install the Linux DEB packages: libssl-dev, libxml2-dev, libcurl4-openssl-dev, lzma, lzma-dev, and liblzma-dev.

Other distributions will have similar packages. The devtools will indicate what Linux installs are missing.

As the user with sudo access and which will run RStudio, set the JAVA_HOME variable. (This is the directory above the bin in, for example, /usr/local/jdk1.8.0_45/bin.)

#Windows Prerequisites

Do a javareconf for R using "R CMD javareconf". Make sure Java 1.8 is registered with R.

Windows installations will also need any extras needed for compiling R packages.

#R Prerequisites

The devtools and rJava package need to be installed.

>install.packages("devtools", dependencies=TRUE)

>install.packages("rJava", dependencies=TRUE)

#Linux Install

>devtools::install_github("GeneSurvey/TCGAGeneReport")

#Windows Install

If you are setup for multiarchitecture compiles in Windows, you can remove the args parameter.
>devtools::install_github("GeneSurvey/TCGAGeneReport", args="--no-multiarch")

#Useful Tests

This can be used to test your rJava install and Java version.

>library(rJava)

>.jinit()

>J("java.lang.System")$getProperty("java.version")
