function ctx_vol = ctx_read_dicomfiles(dcminfo)

[vol,M] = Convert_Dicom_Images(dcminfo);
ctx_vol = mgh2ctx(vol,M);
