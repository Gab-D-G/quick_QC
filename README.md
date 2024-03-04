# Installation

1. Clone this repository, or download all the files and put them in the same folder.
2. Create a python environment with proper dependencies: 

```sh
module load anaconda # need some version of anaconda/miniconda
conda create -n quick_qc python=3.9
conda activate quick_qc
pip install rabies==0.5.1 bruker
```

# Usage

Load dependencies and activate python environment:

```sh
module load anaconda FSL ANTs minc-toolkit-v2 minc-toolkit-extras
conda activate quick_qc
```

### Command line interface --help
```
./quick_qc.py --help
usage: quick_qc.py [-h] [--raw_path RAW_PATH] [--isnii] [--scan_id SCAN_ID] [--reco_id RECO_ID] [--TR TR] [--cutoff CUTOFF] [--ica_dim ICA_DIM]
                   output_folder

Parser to handle testing using token data.

positional arguments:
  output_folder        the output folder path.
                       

options:
  -h, --help           show this help message and exit
  --raw_path RAW_PATH  path to bruker raw folder.
  --isnii              Select this parameter if a nifti file was provided for --raw_path 
                       instead of bruker raw data. No conversion will be run. 
                       (default: False)
                       
  --scan_id SCAN_ID    scan id to select in raw folder.
  --reco_id RECO_ID    reco id if multiple reconstructions. Provide 1 otherwise.
                       (default: 1)
                       
  --TR TR              TR in seconds.(default: 1.0)
                       
  --cutoff CUTOFF      Can specify a cutoff range to select only part of the timeseries 
                       (e.g. '0,120' will select the frames from index 0 to 120, i.e. 2min 
                       of acquisition with a TR of 1.0 second). 'none' takes the full timeseries. 
                       (default: none)
                       
  --ica_dim ICA_DIM    Dimensionality for ICA. Automatic estimation if 0 is provided.
                       (default: 0)
```

### Example 
```sh
quick_qc.py path_to_out_folder/ --raw_data path_to_bruker_raw/ --scan_id 4 --TR 1.0 --cutoff 0,120 --ica_dim 20
```

# Consulting the outputs

* **inho_cor.png**: consult to make sure the brain mask is well registered onto the EPI
* **melodic_masked.ica**: consult the MELODIC html report by clicking report/00index.html
* **DR_fit/**: This folder contains the dual regression outputs. The network fits can be inspected.
