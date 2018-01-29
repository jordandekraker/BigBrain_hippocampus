Scripts for analysis of BigBrain hippocampal blocks (found on graham at /project/6007967/jdekrake/BigBrain/).

Note when using /unfolded.?/indexed/indexes_mapping.mat: this file contains lists of indices idxgm (archicortex) and idxdg (dentate gyrus) voxels, which can be embedded into a 3D logical matrix with size(sz) (e.g. test = zeros(sz); test(idxgm) = 1;). All other files in /unfolded.?/indexed will contain lists of values of the same size as either idxgm or idxdg. Note also that to write data to a nifti file with the same header as the original file, use the variable origheader.img.

TODO:
 - correct isovolume solution for the dentate gyrus
 - get get gyrification index
