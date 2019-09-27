* * * 
## How to verify this pipeline installation on your computer
Use verification testing to assure that the runtime environment and the current version produce the expected output using this repository's data.
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

### 3. Change directory to Feature_Prioritization_Pipeline/test

```
cd Feature_Prioritization_Pipeline/test
```

### 4. Start the verification test from the command line

```
make verification_tests
```

### 5. The output files will be compared with the Feature_Prioritization_Pipeline/data/verification/BENCHMARK_XX... data
* Each Benchmark will report PASS or FAIL and list the names of files producing differences (if any).
* Note that the files generated will be erased after each Benchmark test.

