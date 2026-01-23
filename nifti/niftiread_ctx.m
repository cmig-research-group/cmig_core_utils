function vol_ctx = niftiread_ctx(fname_nii)

[vol M] = niftiread_amd(fname_nii);
vol_ctx = mgh2ctx(vol,M);

