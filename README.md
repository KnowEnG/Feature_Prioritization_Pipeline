# KnowEnG's Feature Prioritization Pipeline
This is the Knowledge Engine for Genomics (KnowEnG), an NIH, BD2K Center of Excellence, Feature Prioritization Pipeline.

This pipeline **ranks** the rows of a given spreadsheet, where spreadsheet's rows correspond to feature-labels and columns correspond to sample-labels. The ranking is based on correlating feature expression data against response data.

There are two prioritization methods:


| **Options**                                        | **Method**                           | **Parameters**            |
| -------------------------------------------------- | -------------------------------------| ------------------------- |
| Simple Correlation                                 | simple correlation                          | correlation               |
| Bootstrap Correlation                              | bootstrap sampling correlation       | bootstrap_correlation     |


Note: all of the correlation methods mentioned above use the Pearson, t-test, or edgeR correlation measure method.

* * * 
## How to run this pipeline with Our data
* * * 

### 1. Clone the Feature_Prioritization_Pipeline Repo
```
 git clone https://github.com/KnowEnG/Feature_Prioritization_Pipeline.git
```

### 2. Install the following (Ubuntu or Linux)
  ```
 apt-get install -y python3-pip
 apt-get install -y libblas-dev liblapack-dev libatlas-base-dev gfortran
 pip3 install numpy==1.11.1
 pip3 install pandas==0.18.1
 pip3 install scipy==0.18.0
 pip3 install scikit-learn==0.17.1
 apt-get install -y libfreetype6-dev libxft-dev
 pip3 install matplotlib==1.4.2
 pip3 install pyyaml
 pip3 install knpackage
 # the following lines are required only for edgeR
 export R_BASE_VERSION=3.6.1
 apt-get update && apt-get install -y apt-transport-https
 echo "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" > \
    /etc/apt/sources.list.d/r.list
 apt-key adv --keyserver keyserver.ubuntu.com \
    --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
 apt-get update && apt-get install -y \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-*
 pip3 install Cython==0.29.13
 pip3 install feather-format==0.3.1
 Rscript Feature_Prioritization_Pipeline/r_src/installation.R
```

### 3. Change directory to Feature_Prioritization_Pipeline

```
cd Feature_Prioritization_Pipeline
```

### 4. Change directory to test

```
cd test
```
 
### 5. Create a local directory "run_dir" and place all the run files in it
```
make env_setup
```

### 6. Use one of the following "make" commands to select and run a prioritization option:


| **Command**                        | **Option**                                        | 
|:---------------------------------- |:------------------------------------------------- | 
| make run_pearson          | pearson correlation                                       |
| make run_bootstrap_pearson | bootstrap sampling with pearson correlation                    |
| make run_t_test          | t-test correlation                                       |
| make run_bootstrap_t_test | bootstrap sampling with t-test correlation                    |
| make run_edgeR          | edgeR correlation                                       |
| make run_bootstrap_edgeR | bootstrap sampling with edgeR correlation                    |

 
* * * 
## How to run this pipeline with Your data
* * * 

__***Follow steps 1-3 above then do the following:***__

### * Create your run directory

 ```
 mkdir run_directory
 ```

### * Change directory to the run_directory

 ```
 cd run_directory
 ```

### * Create your results directory

 ```
 mkdir results_directory
 ```
 
### * Create run_paramters file  (YAML Format)
 ``` 
Look for examples of run_parameters in ./Feature_Prioritization_Pipeline/data/run_files/zTEMPLATE_FP_BENCHMARKS.yml
 ```
### * Modify run_paramters file  (YAML Format)
```
set the spreadsheet and response data file names to point to your data
```

### * Run the Feature Prioritization Pipeline:

  * Update PYTHONPATH enviroment variable
   ``` 
   export PYTHONPATH='../':$PYTHONPATH    
   ```
   
  * Run
   ```
  python3 -m knfeatureprioritization.feature_prioritization -run_directory ./ -run_file zTEMPLATE_FP_BENCHMARKS.yml
   ```

* * * 
## Description of "run_parameters" file
* * * 

| **Key**                    | **Value**                             | **Comments**                                                            |
| -------------------------- | ------------------------------------- | ----------------------------------------------------------------------- |
| method                     | correlation or  bootstrap_correlation | Choose feature prioritization method                                    |
| correlation_measure        | pearson or t_test or edgeR            | Choose correlation measure method                                       |
| spreadsheet_name_full_path | directory+spreadsheet_name            | Path and file name of user supplied feature sets                        |
| phenotype_name_full_path   | directory+response                    | Path and file name of user supplied response file                       |
| results_directory          | directory                             | Directory to save the output files                                      |
| number_of_bootstraps       | 5                                     | Number of random samplings                                              |
| cols_sampling_fraction     | 0.9                                   | Select 90% of spreadsheet columns                                       |
| top_beta_of_sort           | 100                                   | Number of top features selected                                         |
| max_cpu                    | 4                                     | Maximum number of processors to use in the parallel correlation section |

spreadsheet_name = CCLE_Expression_ensembl.df</br>
response_name = CCLE_drug_ec50_cleaned_NAremoved_pearson.txt

* * * 
## Description of Output files saved in results directory
* * * 

* Any method saves separate files per response with name {response}\_{method}\_{correlation_measure}\_{timestamp}\_viz.tsv. Features are sorted in descending order based on `visualization_score`. </br>  

 | **Response**  | **Feature_ID** | **quantitative_sorting_score** | **visualization_score** | **baseline_score** |
 |:-------------:|:--------------:|:------------------------------:|:-----------------------:|:------------------:|
 |   response 1  |   feature 1    |    float                       |    float                |   float            | 
 |    ...        |   ...          |    ...                         |    ...                  |   ...              | 
 |   response 1  |   feature n    |    float                       |    float                |   float            | 


* Any method saves sorted features for each response with name ranked_features_per_response\_{method}\_{correlation_measure}\_{timestamp}\_download.tsv.

 |**Ranking**| **response 1**                  |**response 2**                  |**...**|**response n**                  |
 |:---------:| :------------------------------: |:------------------------------: | :---: |:-------------------------------:|
 |1          | feature </br> (most significant) |feature </br> (most significant) |...    |feature </br> (most significant) |
 |...        |...                               | ...                             |...    |...                              |
 |n          |feature </br> (least significant) |feature </br> (least significant)|...    |feature </br> (least significant)|
 
 
 
* Any method saves spreadsheet with top ranked features per response with name  top_features_per_response\_{method}\_{correlation_measure}\_{timestamp}\_download.tsv.

 |**Features**            | **response 1**      |**...**|**response n**       |
 | :--------------------: |:--------------------:| :---: |:--------------------:|
 | feature 1              |1/0                   |...    |1/0                   |
 | ...                    |...                   |...    |...                   |
 | feature n              | 1/0                  |...    |1/0                   |
