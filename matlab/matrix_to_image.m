
addpath('C:\Program Files\MATLAB\R2021b\spm12'); %based on SPM
x = readmatrix('data\output\all_schaefersc_v4_icc_matrix.csv');
image = 'data\CATD\derivatives\rest_processed\Schaefer2018_400Parcels_7Networks_order_wBascSC_FSLMNI152_2mm.nii';
writetoimage(sum(x), image,452,'out_icc.nii')
