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

## Introductory Example
Make sure you can repeat the example in [EssRegVignette_plainER.pdf](https://github.com/Hanxi-002/EssReg/blob/main/EssRegVignette_plainER.pdf) first. 
Example data: `x` is myeloid and fibroblasts from scRNA-seq from Bob Lafyatis. `y` vector is MRSS skin score. (There are links embedded into the document. Download and open the pdf on your machine to enable interactive links.) This data is linked within the vignette pdfs, but we also provide links here: [x](https://pitt-my.sharepoint.com/:x:/r/personal/xiaoh_pitt_edu/Documents/Research_Files/EssReg/x.csv?d=wcf04e38daa894022b0c2179e4569a141&csf=1&web=1&e=Z0ZcDC), [y](https://pitt-my.sharepoint.com/:x:/r/personal/xiaoh_pitt_edu/Documents/Research_Files/EssReg/y.csv?d=w60482b81ce7e4676b471582c9345ef5a&csf=1&web=1&e=PvgNsP). This example can be run locally. <br>

## Pipeline Example
Once you have reproduced the results in the introductory example, try to run the pipeline using the tutorial found in [EssRegVignette_pipeline.pdf](https://github.com/Hanxi-002/EssReg/blob/main/EssRegVignette_plainER.pdf). (There are links embedded into the document. Download and open the pdf on your machine to enable interactive links.)
The data used here is the same as that in the introductory example. This example can be run locally. <br>

## Detailed Vignette
A comprehensive vignette is available in raw format in the main directory of the package as `EssRegVignette.html`. To view it, go to [this link](https://hanxi-002.github.io/). Because of the size of the data, we do not recommend trying to replicate these results as you will most likely need to submit a job to the cluster. For now, we leave the data in the `data/` folder, but THIS FEATURE MAY BE REMOVED IN THE NEAR FUTURE. In all honestly, being able to reproduce the pipeline example above is sufficient to check your understanding and familiarize yourself with the methods.

## Pipeline Description
### Concept

ER Step 1, coarse grid delta search with grids 0.0001 - 0.001, 0.001 - 0.01, 0.01-0.1, 0.1-1. <br>
ER Step 2, K-CV of the 4 optimal deltas from step 1.  <br>
ER Step 3, find grid delta search from the delta chosen by user in step 2. <br>
ER Step 4, K-CV of the optimal delta from step 3 and multiple lambda values. <br>
ER Step 5, after confirming for the optimal delta and optimal lambda, re-run KCV with larger number of reps. (OPTIONAL)

### Implementation
**Data Format**: Please make sure both x AND y have row names and column names (see vignettes for more information).

All of these functions can be run locally (may be prohibitively slow), in an interactive session onDemand, or using a bash script on the cluster. Below we provide the R function calls, which you can type into RStudio or whatever IDE you are using. If using bash, then you will need to create an R script that contains the function calls below. 

The templates can be found in the `data/` folder.

Summary: <br>
1. For ER Steps 1 & 2, edit `config_template1.yaml`, which is the configuration file for the function `pipelineER1`:
```R
pipelineER1(yaml_path = "config_template1.yaml") ## change path
```
2. For ER Steps 3 & 4, edit `config_template2.yaml`, which is the configuration file for the function `pipelineER2`:
```R
pipelineER2(yaml_path = "config_template2.yaml") ## change path
```
3. For ER Step 5, edit `config_template3.yaml`, which is the configuration file for the function `pipelineER3`:
```R
pipelineER3(yaml_path = "config_template3.yaml") ## change path
```
		
## References
[Essential Regression - a generalizable framework for inferring causal latent factors from multi-omic human datasets](https://www.biorxiv.org/content/10.1101/2021.05.03.442513v2)

## FAQ

**I tried re-installing the package, but I am still getting the same errors I was seeing before the re-installation**
- If you are using RStudio, you will want to restart your R session. This can be done by going to `Session -> Restart R`. Be aware that you will need to reload your data/objects/libraries after restarting. 

**I wanted to try some new parameters so I re-ran the code, but my results are not appearing/not different from the last run**
- First, make sure you are working in the correct directory and that your output pathway is correctly specified. If you are sure that you are looking in the correct directory, then it may be that you already had output created from a prior run. To avoid re-doing time-consuming calculations, Essential Regression performs file existence checks and does not recompute costly steps if they have been previously done. However, changing parameter values does not override these checks. If you want to re-run code fully, you will need to do one of the following: delete the old output directory, rename the old output directory, use a different output pathway.
