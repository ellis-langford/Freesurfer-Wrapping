<div align="center">
  <img src="./assets/freesurfer_logo.png" width="700">
  <br><br>
  <p align="center"><strong>FreeSurfer Wrapping: An open-source software suite for processing human brain MRI</strong></p>
</div>

<div align="center" style="display: flex; justify-content: center; gap: 10px; flex-wrap: wrap; margin-top: 10px;">
  <a href="https://profiles.ucl.ac.uk/101480-ellis-langford"><img src="https://custom-icon-badges.demolab.com/badge/UCL Profile-purple?logo=ucl" alt="UCL Profile"></a>
  <a href="https://orcid.org/0009-0006-1269-2632"><img src="https://img.shields.io/badge/ORCiD-green?logo=orcid&logoColor=white" alt="ORCiD"></a>
  <a href="https://github.com/ellis-langford"><img src="https://img.shields.io/badge/GitHub-%23121011.svg?logo=github&logoColor=white" alt="GitHub"></a>
  <a href="https://uk.linkedin.com/in/ellis-langford-8333441ab"><img src="https://custom-icon-badges.demolab.com/badge/LinkedIn-0A66C2?logo=linkedin-white&logoColor=fff" alt="LinkedIn"></a>
</div>

## Introduction

This pipeline is a wrapping of the FreeSurfer pipeline which performs surface and subcortical segmentation of brain MRI images to generate surface (.stl) files for meshing.


## Requirements

To successfully run the FreeSurfer Wrapping pipeline, please ensure the following requirements are met:

**Ubuntu 22.04 + Docker 27.3.1 + Python 3.10**<br>
*(other versions may be compatible but have not been tested)*


## Installation & Quick Start

To install the necessary components for FreeSurfer Wrapping, please follow the steps below:

- Either, pull the docker image from GitHub container registry:

  ```bash
  docker pull ghcr.io/ellis-langford/freesurfer:v1
  ```

- Or clone the code from the GitHub repo and build image yourself:
  
  ```bash
  git clone https://github.com/ellis-langford/FreeSurfer-Wrapping.git
  cd FreeSurfer-Wrapping
  docker build -t ghcr.io/ellis-langford/freesurfer:v1 .
  ```
  
- Ensure to add your own FreeSurfer license file to /app/license.txt

- Lauch a docker container from the FreeSurfer Wrapping docker image:
  
  ```bash
  docker run -it -v /path/to/data:/path/to/data ghcr.io/ellis-langford/freesurfer:v1 bash
  ```

- Edit the example properties file to suit your requirements
  
  ```bash
  nano example_properties_file.json
  ```

- Navigate to your chosen output directory:
  
  ```bash
  cd /outputdir
  ```

- Run the pipeline:
  
  ```bash
  python3.10 /app/core/freesurfer.py --subject_id subjectID --input_im /path/to/input/image --props_fpath /path/to/properties/file


## Pipeline Options

After running these commands, the FreeSurfer Wrapping pipeline will be run according to the options in the properties file. The pipeline steps available include:

1. `recon-all`: Perform FreeSurfer subcortical segmentation
   > *parallelise:* Run freesurfer processing in parallel. Default is True.<br>
   > *no_randomness:* Process with fixed random seeds to ensure reproducibility. Default is True.<br>
   > *big_vents:* Use if subject has enlarged ventricles from atrophy. Default is False.<br>
   > *large_FOV:* Use if subject has a field of view > 256. Default is False.<br>
   > *optional_flags:* String of comma seperated flags to include in command.<br>
2. `Convert segmentation to NIfTI`: Convert aseg.mgz to aseg.nii.gz.
3. `Convert T1 to NIfTI`: Convert T1.mgz to T1.nii.gz.
4. `Binarise region labels`: Create region speicifc segmentations by selecting label indices
5. `Register atlas labels to subject space`: Register MNI-ICMB152 atlas labels to subject space
6. `Split brainstem labels`: Split brainstem labels into left and right aided by atlas
7. `Generate regional surface files`: Tesselate, smooth and convert regions to .stl

Input options for running the pipeline include:
1. T1 image for full FreeSurfer pipeline processing:
    > *--input_im* Path to input T1.nii.gz
2. A directory containing already processed FreeSurfer outputs to carry out post processing:
    > *--freesurfer_outputs* Path to FreeSurfer output directory
3. A comma seperated list of paths to segmentation files for post-processing:<br>
   Note: if both wholebrain and ventricle segmentations are provided, a global segmentation and surface are generated
    > *--segmentations* A comma seperated list of paths to segmentation files

## Data Preparation

No additional data preparation is required.


## Output Structure

The output directory structure is as follows:

```
Output directory
├── inputs
├── logs
├── outputs
├── fs_outputs
├── results.txt
└── errors.txt
```
- `inputs:` contains a copy of the input images
- `logs:` contains a plugin log (log.txt) and a record of the inputs and parameters (options.txt)
- `outputs:` contains the final output images
- `fs_outputs`: contains the raw output folder produced by the FreeSurfer pipeline
- `results.txt:` only produced if the pipeline executes successfully
- `errors.txt:` only produced if the pipeline fails to execute successfully (contains error info)


## Citation

```
@ARTICLE{xxxxxxxx,
  author={Langford E},
  journal={}, 
  title={}, 
  year={},
  volume={},
  number={},
  pages={},
  doi={}}
```

## Useful links

[![](https://img.shields.io/badge/Software-FreeSurfer-orange)](https://surfer.nmr.mgh.harvard.edu/fswiki)


## References

1. *B. Fischl, “FreeSurfer,” NeuroImage, vol. 62, no. 2, pp. 774–781, 2012.*
