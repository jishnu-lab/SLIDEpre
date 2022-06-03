library(EssReg)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))
if(is.na(cores)) cores <- detectCores()
# cores <- 6
registerDoParallel(cores)
cat('number of cores using', cores, '. . .\n')

yaml_path = '/ihome/djishnu/xiaoh/Lafyatis/All_Cell_Type/HER_041422/ER_Config_Step3_4.yml'
pipelineER2(yaml_path)

