# Installation

1. Clone this repository, or download all the files and put them in the same folder.
2. Create a python environment with proper dependencies: 

```sh
module load anaconda # need some version of anaconda/miniconda
conda create -n quick_qc python=3.9
conda activate quick_qc
pip install rabies==0.5.1
```

# Usage

Load dependencies and activate python environment:
'''
module load anaconda FSL ANTs minc-toolkit-v2 minc-toolkit-extras
conda activate quick_qc
'''

Run the script:

```sh
bash quick_qc.sh $bold_file $out_folder $TR $cutoff $ica_dim
```

* bold_file: the path to the fMRI nifti file
* out_folder: output folder path
* TR: specify the TR in seconds (e.g. 1.0)
* cutoff: can specify a cutoff range to select only part of the timeseries (e.g. '0,120' will select the frames from index 0 to 120, i.e. 2min of acquisition with a TR of 1.0 second). To take the full timeseries, simply write 'none'
* ica_dim: provide the dimensionality for the ICA decomposition (# of components to derive). Provide 0 to conduct MELODIC's automatic dimensionality estimation.

Example:
```sh
bash quick_qc.sh bold_file.nii out_folder 1.0 '0,120' 0
```

# Consulting the outputs

* **inho_cor.png**: consult to make sure the brain mask is well registered onto the EPI
* **melodic_masked.ica**: consult the MELODIC html report by clicking report/00index.html
* **DR_fit/**: This folder contains the dual regression outputs. The network fits can be inspected.
