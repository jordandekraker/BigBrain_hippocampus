% Upsample back to 40um (from 80um versions) using linear interpolation and
% NaN edge correction. Applied to the labelmap, and each gradient.
clear all; close all; 

for LR = 'lr'

%% the labelmap

mkdir(sprintf('unfolded.%s',LR));

origheader = load_untouch_nii(sprintf('hippocampus.%s.40um.nii.gz',LR));
origheader.img = [];
ds_lbl = load_untouch_nii(sprintf('downsampled/Unfolded.%s/labelmap.%s.nii.gz',LR,LR));
load(sprintf('downsampled/Unfolded.%s/data.mat',LR));

ds_idxgm = find(ds_lbl.img==2);
ds_idxdg = find(ds_lbl.img==1);
ds_sz = size(ds_lbl.img);

lbl = origheader;
lbl.img = imresize3(ds_lbl.img,2,'nearest');
save_untouch_nii(lbl,sprintf('unfolded.%s/labelmap.nii.gz',LR))
idxgm = find(lbl.img==2);
idxdg = find(lbl.img==1);
sz = size(lbl.img);
clear ds_lbl lbl

mkdir(sprintf('unfolded.%s/indexed',LR));
mkdir(sprintf('unfolded.%s/nifti(discretized)',LR));

save(sprintf('unfolded.%s/indexed/indexes_mapping.mat',LR),'sz', 'ds_sz', 'idxgm', 'idxdg', 'ds_idxgm', 'ds_idxdg', 'origheader');

%% now each of the Laplace gradients

hl=ones(3,3,3);
hl = hl/26; hl(2,2,2) = 0;
nanedge_correction_gm = zeros(sz); 
nanedge_correction_gm(idxgm) = 1;
nanedge_correction_gm = imfilter(nanedge_correction_gm,hl,'replicate','conv');

%PD gm
origheader.img = zeros(ds_sz);
origheader.img(ds_idxgm) = Laplace_PD;
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_gm; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
PD_laplace = origheader.img(idxgm);
save(sprintf('unfolded.%s/indexed/PD_laplace.mat',LR),'PD_laplace');
origheader.img = ceil(origheader.img*100); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/PD_laplace.nii.gz',LR));

%AP gm
origheader.img = zeros(ds_sz);
origheader.img(ds_idxgm) = Laplace_AP(1:length(ds_idxgm));
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_gm; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
AP_laplace = origheader.img(idxgm);
save(sprintf('unfolded.%s/indexed/AP_laplace.mat',LR),'AP_laplace');
origheader.img = ceil(origheader.img*100); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/AP_laplace.nii.gz',LR));

%IO laplace gm
origheader.img = zeros(ds_sz);
origheader.img(ds_idxgm) = Laplace_IO;
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_gm; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
IO_laplace = origheader.img(idxgm);
save(sprintf('unfolded.%s/indexed/IO_laplace.mat',LR),'IO_laplace');
origheader.img = ceil(origheader.img*100); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/IO_laplace.nii.gz',LR));

%IO isovolume gm
isovolume = load_untouch_nii(sprintf('downsampled/Unfolded.%s/lamination.%s_layering_depth.nii.gz',LR,LR));
origheader.img = zeros(ds_sz);
origheader.img(ds_idxgm) = isovolume.img(ds_idxgm); 
clear isovolume;
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_gm; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
IO_isovolume = origheader.img(idxgm);
save(sprintf('unfolded.%s/indexed/IO_isovolume.mat',LR),'IO_isovolume');
origheader.img = ceil(origheader.img*10); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/IO_isovolume.nii.gz',LR));

clear nanedge_correction_gm
nanedge_correction_dg = zeros(sz); 
nanedge_correction_dg(idxdg) = 1;
nanedge_correction_dg = imfilter(nanedge_correction_dg,hl,'replicate','conv');

%AP dg
origheader.img = zeros(ds_sz);
origheader.img(ds_idxdg) = Laplace_AP(length(ds_idxgm)+1:end);
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_dg; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
AP_laplace_dg = origheader.img(idxdg);
save(sprintf('unfolded.%s/indexed/AP_laplace_dg.mat',LR),'AP_laplace_dg');
origheader.img = ceil(origheader.img*100); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/AP_laplace_dg.nii.gz',LR));


%IO dg
isovolume = load_untouch_nii(sprintf('downsampled/Unfolded.%s/laminationDG.%s_layering_depth.nii.gz',LR,LR));
origheader.img = zeros(ds_sz);
origheader.img(ds_idxdg) = isovolume.img(ds_idxdg); 
clear isovolume;
origheader.img = imresize3(origheader.img,2,'nearest');
origheader.img = imfilter(origheader.img,hl,'replicate','conv'); %apply averaging filter
origheader.img = origheader.img./nanedge_correction_dg; %this and the previous line amount to linear interpolation, but without including values outside of idxgm (i.e. considers them NaN)
IO_isovolume_dg = origheader.img(idxdg);
save(sprintf('unfolded.%s/indexed/IO_isovolume_dg.mat',LR),'IO_isovolume_dg');
origheader.img = ceil(origheader.img*10); %for visualization
save_untouch_nii(origheader,sprintf('unfolded.%s/nifti(discretized)/IO_isovolume_dg.nii.gz',LR));

clearvars -except LR
end
clear
