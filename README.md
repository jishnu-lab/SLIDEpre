# EssReg
This repo is currently for in lab personnel only. 
## Installation
```R
devtools::install_github("Hanxi-002/EssReg"
                         ,ref="main"
                         ,auth_token = "tokenstring"
                         )
```

[How to generate an authentification token?](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token)
## Pipeline Description
### Concept

ER Step 1, coarse grid delta search with grids 0.0001 - 0.001, 0.001 - 0.01, 0.01-0.1, 0.1-1. <br>
ER Step 2, K-CV of the 4 optimal deltas from step 1.  <br>
ER Step 3, find grid delta search from the delta chosen by user in step 2. <br>
ER Step 4, K-CV of the optimal delta from step 3 and multiple lambda values. <br>
ER Step 5, after confirming for the optimal delta and optimal lambda, re-run KCV with larger number of reps. 

	
### Implementation
**Data Format**: Please make sure both x AND y have row names and column names!! (see exmaple or bench mark folder for examples of x and y csv files.)

TLDR: <br>
1. First create a folder to contain the output and use the path as outpath in the config files. <br>
2. Edit `ER_Config_Step1_2.yml`. <br>
`runER_shell.sh` (needs edits) --> `run_ERPipeS1.R` (needs edits) --> `pipelineER1.R` --> read box plot to choose best delta. <br>
3. Edit `ER_Config_Step3_4.yml`. <br>
4. `runER_shell.sh` (needs edits) --> `run_ERPipeS2.R`(needs edits) --> `pipelineER2.R` --> read box plot to choose best lambda. <br>
5. Edit `ER_Config_Step5.yml`. <br>
6. `runER_shell.sh` (needs edits) --> `run_ERPipeS3.R`(needs edits) --> `pipelineER3.R`. <br>

Details: <br>
1. For ER Step 1 & 2, run function `pipelineER1`. See `run_ERPipeS1.R` in exmaple folder to call `pipelineER1`. <br>
	Input: `ER_Config_Step1_2.yml` <br>
	Output: 4 folders, one corresponds to each delta grid, with correlation heat maps. `er_s1_output.RDS` which contains all 4 ER output from the first 4 delta grid search. `delta_selection_boxplot.pdf` that shows the cross validation results with permutation. 
2. For ER Step 2 & 3, fun function `pipelineER2`. See `run_ERPipeS2.R` in example folder to call `pipelineER2`. <br>
	Input: `ER_Config_Step3_4.yml` <br>
	Output:`er_s1_output.RDS` which contains the ER output from step 3. `lambda_selection_boxplot.pdf` that shows the cross validation results with permutation. 
	
		
## References
[Essential Regression - a generalizable framework for inferring causal latent factors from multi-omic human datasets](https://www.biorxiv.org/content/10.1101/2021.05.03.442513v2)
## ER Example
Make sure you can repeat the toy example first. A separate Readme file is included in this folder. <br>
Toy example data: myeloid and fibroblasts from scRNA-seq from Bob Lafyatis. Y vector is MRSS skin score. <br>
## Benchmark plainER
plainER is identifiable which means determinstic parameters should produce the exact same resutls. This is an insanity check. 
