###INPUTS

# example: bash quick_qc.sh bold_file.nii out_folder 1.0 '0,120'

bold_file=$1
output_folder=$2
TR=$3
cutoff=$4
ica_dim=$5

masking_thresh=1.2
DR_ON=TRUE
###

# get path of the script
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

template_anat=$parent_path/EPI_template.nii.gz
template_mask=$parent_path/EPI_brain_mask.nii.gz
SM_IC_file=$parent_path/somatomotor_IC.nii.gz
DMN_IC_file=$parent_path/DMN_IC.nii.gz


echo "Process timeseries"
mkdir -p $output_folder
python $parent_path/process.py $bold_file $TR $cutoff $output_folder # takes ~6 seconds to run


echo "Bias correction"
bias_cor_out=$output_folder/bias_cor
mkdir -p $bias_cor_out
input=$output_folder/median.nii.gz
ImageMath 3 ${bias_cor_out}/null_mask.nii.gz ThresholdAtMean $input 0
ImageMath 3 ${bias_cor_out}/thresh_mask.nii.gz ThresholdAtMean $input $masking_thresh

N4BiasFieldCorrection -d 3 -s 4 -i $input -b [20] -c [200x200x200,0.0] -w ${bias_cor_out}/thresh_mask.nii.gz \
-x ${bias_cor_out}/null_mask.nii.gz -o ${bias_cor_out}/N4.nii.gz
DenoiseImage -d 3 -i ${bias_cor_out}/N4.nii.gz -o ${bias_cor_out}/denoise.nii.gz


echo "Registration"

reg_out=$output_folder/reg_out
mkdir -p $reg_out
reg_type="--linear-type rigid --skip-nonlinear"
antsRegistration_affine_SyN.sh $reg_type --convergence 1e-7 ${bias_cor_out}/denoise.nii.gz ${template_anat} ${reg_out}/tomodel

antsApplyTransforms -d 3 -i ${template_mask} -t [ ${reg_out}/tomodel0GenericAffine.mat,1 ] \
-o ${reg_out}/mask_resampled.nii.gz -n GenericLabel -r ${input} --verbose


mkdir -p $output_folder/resampled_ICs
antsApplyTransforms -d 3 -i $SM_IC_file -t [ ${reg_out}/tomodel0GenericAffine.mat,1 ] \
-o $output_folder/resampled_ICs/SM_IC.nii.gz -n Linear -r ${input} --verbose

antsApplyTransforms -d 3 -i $DMN_IC_file -t [ ${reg_out}/tomodel0GenericAffine.mat,1 ] \
-o $output_folder/resampled_ICs/DMN_IC.nii.gz -n Linear -r ${input} --verbose

echo "Dual regression"
python $parent_path/process2.py $output_folder $DR_ON

echo "MELODIC"

# way faster to compute with the mask, ~12sec
melodic -i $output_folder/processed_img.nii.gz -o $output_folder/melodic_masked.ica \
--report -v -m ${reg_out}/mask_resampled.nii.gz --tr=$TR --bgimage=$output_folder/median.nii.gz --seed=1 -d $ica_dim

# generate melodic at the end within mask in case registration failed
melodic -i $output_folder/processed_img.nii.gz -o $output_folder/melodic_nomask.ica \
--report -v --nomask --tr=$TR --bgimage=$output_folder/median.nii.gz --seed=1 -d $ica_dim
