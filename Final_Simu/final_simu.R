#!/usr/bin/R

suppressMessages( library(magrittr) )
suppressMessages( library(argparser) )

p <- arg_parser("Simulation on Group lasso compared with multiple other methods") %>%
  # add_argument(
  #   arg = "--seed",
  #   help = "Random seed",
  #   type = "integer"
  # ) %>%
  # add_argument(
  #   arg = "--samples",
  #   help = "Number of samples",
  #   default = 100,
  #   short = "-N"
  # ) %>%
  # add_argument(
  #   arg = "--baseline_dim",
  #   help = "Number of baseline covariates",
  #   default = 2,
  #   short = "-m_X"
  # ) %>%
  # add_argument(
  #   arg = "--treatment_dim",
  #   help = "Number of treatment covariates",
  #   default = 2,
  #   short = "-m_W"
  # ) %>%
  # add_argument(
  #   arg = "--gene_dim",
  #   help = "Number of gene covariates",
  #   default = 10,
  #   short = "-m_G"
  # ) %>%
  # add_argument(
  #   arg = "--out",
  #   help = "Output file",
  #   default = ""
  # ) %>%
  # add_argument(
  #   arg = "--main_nonzero",
  #   help = "Proportion of nonzero among all nonzero main effects",
  #   default = .1,
  #   short = "-MN"
  # ) %>%
  # add_argument(
  #   arg = "--interaction_nonzero",
  #   help = "Proportion of nonzero among all nonzero interaction effects",
  #   default = .1,
  #   short = "-IN"
  # ) %>%
  # add_argument(
  #   arg = "--both_nonzero",
  #   help = "Proportion of both nonzero among all nonzero interaction effects when hierachical relationship doesn't exist",
  #   default = .1,
  #   short = "-BN"
  # ) %>%
  add_argument(
    arg = "--SNR",
    help = "Signal Noise Ratio",
    type = "integer",
    default = 10
  ) %>%
  add_argument(
    arg = "--Cov_type",
    help = "Data type of gene covariates",
    default = "SNP",
    short="-type"
  ) %>%
  add_argument(
    arg = "--Hierarchy",
    help = "Whether hierarchical relationship exist",
    default = TRUE,
    type = "logical",
    short="-H"
  ) 
  argv <- parse_args(p)

if (length(argv) < 3) {
  stop("Must have at least 12 covariates")
}


#working directory, need to be changed
dir<-"/home/fas/zhao/wd262/project/predictive" 

#need to import library(dr) and library(MASS)
source(sprintf("%s/cv3.R",dir))
source(sprintf("%s/func3.R",dir))
source(sprintf("%s/opt3.R",dir))
source(sprintf("%s/sim4.R",dir))

library(rpart)
library(randomForest)
library(BMA)
library(leaps)
#library(BoomSpikeSlab)
library(SIS)

print(argv)