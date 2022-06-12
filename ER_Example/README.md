# ER Toy Example

Please carefully review the vignette before trying this example. Some of the paths will needs to be changed for your system and file locations.

## Input Data
1. Input X data: `x_prescale.csv`. <br>
	Input Y data: `check_y.csv`.
	
## Scripts and Results Description 

### Step1 & Step2 (Scripts):

1. The *pipelineER1* function runs the first **2** steps of the ER Pipeline. To call the pipelineER1 function, please see the `run_ERPipeS1.R`.
2. *pipelineER1* function takes in a string that is a path to a yml file. See `pipeline1.yaml` for an example of the configuration file. 

### Step1 & Step2 (Output Files):
1. All the output files for Step1 and Step2 is in the results_1 folder. The vignette contains information and explanation of what these output files entail. 

### Step3 & Step4 (Scripts):

1. The *pipelineER2* function runs the last **2** steps of the ER Pipeline. To call the pipelineER1 function, please see the `run_ERPipeS2.R`.
2. *pipelineER2* function takes in a string that is a path to a yml file. See `pipeline2.yaml` for an example of the configuration file. 