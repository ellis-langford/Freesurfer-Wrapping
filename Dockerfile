# Base image
FROM ubuntu:22.04

# Define some environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/London

# Set working directory
WORKDIR /app

# Copy codebase into image
COPY ./core /app/core
COPY ./atlas /app/atlas
COPY ./requirements.txt /app/requirements.txt
COPY ./license.txt /app/license.txt
COPY ./README.md /app/README.md

# Update the package list and install required dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    build-essential \
    curl \
    git \
    && apt-get clean

# Download utility functions
RUN git clone https://github.com/ellis-langford/ImageUtils.git utils

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install JupyterLab and libraries
RUN pip install jupyterlab

# ANTsPy setup
RUN pip install antspyx
ENV ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
ENV ANTS_RANDOM_SEED=1

# Download FreeSurfer
RUN cd /usr/local \
    && wget https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-ubuntu22_amd64-7.4.1.tar.gz \
    && tar -xzvf freesurfer-linux-ubuntu22_amd64-7.4.1.tar.gz \
    && rm -rf freesurfer-linux-ubuntu22_amd64-7.4.1.tar.gz

# Initialise FreeSurfer and copy license
ENV FREESURFER_HOME=/usr/local/freesurfer
RUN printf "\n# Freesurfer\nsource $FREESURFER_HOME/SetUpFreeSurfer.sh\n" >> /root/.bashrc
COPY ./license.txt $FREESURFER_HOME

# Add more custom folders to PYTHONPATH
ENV PYTHONPATH="/app/core:/app/utils"

# Expose JupyterLab's default port (8888)
EXPOSE 8888

# Create a working directory
WORKDIR /

# Set the command to run JupyterLab when the container starts
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]