function ctx_save_nifti(ctx_struct, fname)
% ctx_save_nifti - Write a NIFTI file from a CMIG ctx structure
%
% Input Arguments
%   ctx_struct: ctx structure with image data to be written to NIFTI file
%   fname: Name of NIFTI file to be saved

info.Version = 'NIfTI1';
info.Description = 'NIFTI from CMIG ctx';  
info.SpaceUnits = 'Millimeter';
info.TimeUnits = 'Second';

info.Datatype = 'single';
ctx_struct.imgs = single(ctx_struct.imgs);

[rows, cols, slices, frames] = size(ctx_struct.imgs);
info.ImageSize = [rows, cols, slices, frames];
info.SliceCode = 'Unknown';

info.FrequencyDimension = 0;
info.PhaseDimension = 0;
info.SpatialDimension = 0;

% Enforce positive Qfactor (slice direction)
info.Qfactor = 1;
M_RAS = M_LPH_TO_RAS * ctx_struct.Mvxl2lph;
slice_dir = M_RAS(1:3,3) / norm(M_RAS(:,3));
if slice_dir(3) < 0
  M_flip = M_RAS;
  slice_thk = norm(M_RAS(:,3));
  M_flip(1:3,4) = M_RAS(1:3,4) + (slices * slice_dir * slice_thk);
  M_flip(1:3,3) = -M_flip(1:3,3);
  M_RAS = M_flip;
  ctx_struct.imgs = flip(ctx_struct.imgs, 3);
end

info.TransformName = 'Sform';
info.Transform.Dimensionality = 3;

info.Transform.T = M_RAS';

pixdim_r = norm(M_RAS(:,1));
pixdim_c = norm(M_RAS(:,2));
pixdim_s = norm(M_RAS(:,3));
info.PixelDimensions = [pixdim_r pixdim_c pixdim_s 1];

gz_flag = ~isempty(regexpi(fname, 'gz'));
[filepath, name, ext] = fileparts(fname);
if gz_flag
  [~, name2, ext] = fileparts(name);
  fname = fullfile(filepath, name2);
  niftiwrite(single(ctx_struct.imgs), fname, info, 'Compressed', true);
else
  fname = fullfile(filepath, name);
  niftiwrite(single(ctx_struct.imgs), fname, info);
end

end
