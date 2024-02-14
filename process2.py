'''
Plot mask
'''

import sys
import os
import SimpleITK as sitk
import numpy as np

output_folder = os.path.abspath(sys.argv[1])
DR_ON = sys.argv[2]

raw_img = f'{output_folder}/median.nii.gz'
warped_mask=f'{output_folder}/reg_out/mask_resampled.nii.gz'
init_denoise=f'{output_folder}/bias_cor/denoise.nii.gz'

def inho_cor_diagnosis(raw_img,init_denoise,warped_mask,out_dir,figure_format):
    import os
    import pathlib
    import SimpleITK as sitk
    # set default threader to platform to avoid freezing with MultiProc https://github.com/SimpleITK/SimpleITK/issues/1239
    sitk.ProcessObject_SetGlobalDefaultThreader('Platform')
    os.makedirs(out_dir, exist_ok=True)

    import matplotlib.pyplot as plt
    from rabies.visualization import plot_3d, otsu_scaling
    fig,axes = plt.subplots(nrows=3, ncols=3, figsize=(12*3,2*3))

    scaled = otsu_scaling(raw_img)
    axes[0,0].set_title('Raw Image', fontsize=30, color='white')
    #add_filenames(axes[-1,0], {'File':raw_img})
    plot_3d(axes[:,0],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    axes[0,2].set_title('Resampled Mask', fontsize=30, color='white')
    #add_filenames(axes[-1,2], {'Mask File':warped_mask,'EPI File':raw_img})
    plot_3d(axes[:,2],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')
    sitk_mask = sitk.ReadImage(warped_mask,sitk.sitkFloat32)
    # resample mask to match template
    sitk_mask = sitk.Resample(sitk_mask, scaled)
    plot_3d(axes[:,2],sitk_mask,fig=fig,vmin=-1,vmax=1,cmap='bwr', alpha=0.3, cbar=False)

    scaled = otsu_scaling(init_denoise)
    axes[0,1].set_title('Intensity correction', fontsize=30, color='white')
    #add_filenames(axes[-1,1], {'File':init_denoise})
    plot_3d(axes[:,1],scaled,fig=fig,vmin=0,vmax=1,cmap='viridis')

    plt.tight_layout()
    fig.savefig(f'{out_dir}/inho_cor.{figure_format}', bbox_inches='tight')

inho_cor_diagnosis(raw_img,init_denoise,warped_mask,out_dir=f'{output_folder}/',figure_format='png')



'''
Compute DR
'''
if DR_ON=='TRUE':
    from rabies.analysis_pkg.analysis_math import dual_regression
    import matplotlib.pyplot as plt
    import nilearn.plotting
    from rabies.analysis_pkg.diagnosis_pkg.analysis_QC import masked_plot, percent_threshold
    from rabies.visualization import otsu_scaling, plot_3d

    def plot_DR(sitk_img,map,scaled, mask_file, figpath):
        
        # select vmax at 98th percentile value
        vector = map.copy()
        vector.sort()
        vmax = map.max()
        planes=('sagittal', 'coronal', 'horizontal')
        threshold=vmax*0.01
        
        fig,axes = plt.subplots(nrows=3, ncols=1, figsize=(20,2*3))
        plot_3d(axes,scaled,fig,vmin=0,vmax=1,cmap='gray', alpha=1, cbar=False, num_slices=8, planes=planes)
        cbar_list = plot_3d(axes,sitk_img,fig,vmin=-vmax,vmax=vmax,cmap='cold_hot', alpha=1, cbar=True, threshold=threshold, num_slices=8, planes=planes)
        
        for cbar in cbar_list:
            cbar.ax.get_yaxis().labelpad = 35
            cbar.set_label("Beta \nCoefficient", fontsize=12, rotation=270, color='white')
            cbar.ax.tick_params(labelsize=12)
        
        plt.tight_layout()
        fig.savefig(figpath, bbox_inches='tight')

    DR_out=output_folder+'/DR_fit/'
    os.makedirs(DR_out, exist_ok=True)

    SM_file=f'{output_folder}/resampled_ICs/SM_IC.nii.gz'
    DMN_file=f'{output_folder}/resampled_ICs/DMN_IC.nii.gz'
    mask_file=f'{output_folder}/reg_out/mask_resampled.nii.gz'
    bold_file = f'{output_folder}/processed_img.nii.gz'
    median_file = f'{output_folder}/median.nii.gz'

    brain_mask = sitk.GetArrayFromImage(sitk.ReadImage(mask_file, sitk.sitkFloat32))
    volume_indices = brain_mask.astype(bool)

    data_img = sitk.ReadImage(bold_file, sitk.sitkFloat32)
    data_array = sitk.GetArrayFromImage(data_img)
    num_volumes = data_array.shape[0]
    timeseries = np.zeros([num_volumes, volume_indices.sum()])
    for i in range(num_volumes):
        timeseries[i, :] = (data_array[i, :, :, :])[volume_indices]


    for name,IC_file in zip(['SM','DMN'],[SM_file, DMN_file]):
        IC_arr = sitk.GetArrayFromImage(sitk.ReadImage(IC_file, sitk.sitkFloat32))[volume_indices]
        
        DR = dual_regression(IC_arr.reshape(1,-1), timeseries)
        from rabies.utils import recover_3D
        IC_img = recover_3D(mask_file, DR['C'][:,0])
        sitk.WriteImage(IC_img, f'{DR_out}/{name}_fit.nii.gz')
        
        '''
        Plot DR
        '''
        scaled = otsu_scaling(median_file)
        map = DR['C'][:,0].flatten()
        sitk_img = recover_3D(
            mask_file, map)
        
        plot_DR(sitk_img, map,scaled, mask_file, figpath=f'{DR_out}/{name}_plot_unthresholded.png')
        
        threshold = percent_threshold(map)
        mask=np.abs(map)>=threshold # taking absolute values to include negative weights
        mask_img = recover_3D(mask_file,mask)
        masked = sitk_img*mask_img
        
        plot_DR(masked, map,scaled, mask_file, figpath=f'{DR_out}/{name}_plot_thresholded.png')
