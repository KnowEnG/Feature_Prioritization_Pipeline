FROM knowengdev/base_image:08_01_2018

ENV SRC_LOC /home

# Install the latest knpackage
RUN pip3 install -I knpackage redis

# Install edgeR
ENV R_BASE_VERSION 3.6.1
RUN apt-get update && apt-get install -y apt-transport-https
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu trusty-cran35/" > \
    /etc/apt/sources.list.d/r.list
RUN apt-key adv --keyserver keyserver.ubuntu.com \
    --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update && apt-get install -y \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-*
RUN pip3 install Cython==0.29.13
RUN pip3 install feather-format==0.3.1

# Clone from github
RUN git clone https://github.com/KnowEnG/Feature_Prioritization_Pipeline.git ${SRC_LOC} 

# Set up working directory
WORKDIR ${SRC_LOC}

RUN Rscript r_src/installation.R