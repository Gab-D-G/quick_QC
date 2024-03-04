import argparse
def get_parser():
    """Build parser object"""
    parser = argparse.ArgumentParser(
        description=
            "Parser to handle testing using token data.",
        formatter_class=argparse.RawTextHelpFormatter)
    preprocess.add_argument(
        'output_folder', action='store', type=Path,
        help=
            "the output folder path.\n"
            "\n"
        )
    parser.add_argument(
        "--raw_path", dest='raw_path', type=Path,
        help=
            "path to bruker raw folder."
        )
    parser.add_argument(
        "--scan_id", dest='scan_id', type=int,
        help=
            "scan id to select in raw folder."
        )
    parser.add_argument(
        "--reco_id", dest='scan_id', type=int,
        default=1,
        help=
            "reco id if multiple reconstructions. Provide 1 otherwise.\n"
            "(default: %(default)s)\n"
            "\n"
       )
    parser.add_argument(
        "--TR", dest='TR', type=float,
        help=
            "TR in seconds."
        )
    parser.add_argument(
        "--cutoff", dest='cutoff', type=str,
        default='none',
        help=
            ""
        )
    parser.add_argument(
        "--ica_dim", dest='ica_dim', type=str,
        default='0',
        help=
            "Dimensionality for ICA. Automatic estimation if 0 is provided.\n"
            "(default: %(default)s)\n"
            "\n"
       )

    return parser

parser = get_parser()
opts = parser.parse_args()

###INPUTS
output_folder=opts.output_folder
raw_path = opts.raw_path
scan_id=opts.scan_id
reco_id=opts.reco_id
TR=opts.TR
TR = float(TR)
cutoff=opts.cutoff
ica_dim=opts.ica_dim
DR_ON=True
###INPUTS

from rabies.utils import run_command
from utils import *
import os
script_path = os.path.realpath(__file__)
quick_qc_folder = os.path.dirname(script_path)
os.makedirs(output_folder, exist_ok=True)

print("brkraw conversion")
import brkraw as br
rawdata = br.load(raw_path)
rawdata.save_as(scan_id, reco_id, 'brkraw_out', dir=output_folder, ext='nii.gz')
bold_file = f'{output_folder}/brkraw_out.nii.gz'

print("Process timeseries")
preproc_bold(bold_file, TR, cutoff, output_folder)

run_command(f'bash {quick_qc_folder}/resampling.sh {output_folder}', verbose=True)

raw_img = f'{output_folder}/median.nii.gz'
warped_mask=f'{output_folder}/reg_out/mask_resampled.nii.gz'
init_denoise=f'{output_folder}/bias_cor/denoise.nii.gz'
inho_cor_diagnosis(raw_img,init_denoise,warped_mask,out_dir=f'{output_folder}/',figure_format='png')

if DR_ON:
    print("Dual regression")
    dual_regression_fit(output_folder)

print('MELODIC')
command=f'melodic -i {output_folder}/processed_img.nii.gz -o {output_folder}/melodic_masked.ica \
--report -v -m {output_folder}/reg_out/mask_resampled.nii.gz --tr={str(TR)} --bgimage={output_folder}/median.nii.gz --seed=1 -d {ica_dim}'
run_command(command, verbose=True)


# generate melodic at the end within mask in case registration failed
command=f'melodic -i {output_folder}/processed_img.nii.gz -o {output_folder}/melodic_nomask.ica \
--report -v --nomask --tr={str(TR)} --bgimage={output_folder}/median.nii.gz --seed=1 -d {ica_dim}'
run_command(command, verbose=True)

