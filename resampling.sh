###INPUTS

output_folder=$1
input=$output_folder/median.nii.gz
masking_thresh=1.2
###

# get path of the script
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

template_anat=$parent_path/EPI_template.nii.gz
template_mask=$parent_path/EPI_brain_mask.nii.gz
SM_IC_file=$parent_path/somatomotor_IC.nii.gz
DMN_IC_file=$parent_path/DMN_IC.nii.gz

echo "Bias correction"
bias_cor_out=$output_folder/bias_cor
mkdir -p $bias_cor_out
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
