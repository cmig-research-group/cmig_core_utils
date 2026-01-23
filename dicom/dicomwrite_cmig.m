function dicomwrite_cmig(ctx_struct, output_dir, dcm_hdr_struct)
% dicomwrite_cmig - Write DICOM files for a CMIG ctx structure
%
% Input Arguments
%   ctx_struct: ctx structure with image data to be written to DICOM files
%   output_dir: Path to directory where DICOM files will be saved
%   dcm_hdr_struct: MATLAB structure containing DICOM header info

mkdir(output_dir);

mosaic_flag = 0;
if isfield(dcm_hdr_struct, 'ImageType')
  mosaic_flag = ~isempty(regexpi(dcm_hdr_struct.ImageType, 'MOSAIC'));
end

rtdose_flag = 0;
if isfield(dcm_hdr_struct, 'Modality')
  rtdose_flag = ~isempty(regexpi(dcm_hdr_struct.Modality, 'RTDOSE'));
end

multiframe_flag = 0; % Regular (not-enhanced) multiframe
if isfield(dcm_hdr_struct, 'NumberOfFrames')
  if dcm_hdr_struct.NumberOfFrames > 1
    multiframe_flag = 1;
  end
end

enhanced_flag = 0; % Enhanced DICOM (also multiframe)
if isfield(dcm_hdr_struct, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end

[rows, cols, slices, frames] = size(ctx_struct.imgs);
maxvol = max(ctx_struct.imgs(:));
dat = ctx_struct.imgs;
[hc, hv] = hist(dat(find(dat>0)),1000);
cdf_hc = cumsum(hc)/sum(hc);
val1 = min(hv(find(cdf_hc>0.99)));
val2 = min(hv(find(cdf_hc>0.10)));
val3 = min(hv(find(cdf_hc>0.90)));

fprintf('%s -- %s | Writing DICOM images...\n', datestr(now), mfilename);

total_files = size(ctx_struct.imgs,3) * size(ctx_struct.imgs,4);
if multiframe_flag || enhanced_flag
   total_files = 1;
end

rgb_flag = 0;
if isfield(dcm_hdr_struct, 'PhotometricInterpretation')
  if strcmp(dcm_hdr_struct.PhotometricInterpretation, 'RGB')
    total_files = size(ctx_struct.imgs,3);
    rgb_flag = 1;
  end
end

SeriesInstanceUID = dicomuid;

slice_counter = 1;
vol_counter = 1;
for i = 1:total_files

  if slice_counter > size(ctx_struct.imgs,3)
    slice_counter = 1;
    vol_counter = vol_counter + 1;
  end

  metadata = dcm_hdr_struct;
  metadata.SeriesInstanceUID = SeriesInstanceUID;

  %% Image geometry
  Mvxl2lph = ctx_struct.Mvxl2lph;
  dcr = Mvxl2lph(1:3,1)/norm(Mvxl2lph(1:3,1));  
  dcc = Mvxl2lph(1:3,2)/norm(Mvxl2lph(1:3,2));  
  dcs = Mvxl2lph(1:3,3)/norm(Mvxl2lph(1:3,3));  
  slcthk = norm(Mvxl2lph(1:3,3));
  PixelSpacing = colvec(sqrt(sum(Mvxl2lph(:,[1 2]).^2)));
  ImageOrientationPatient = [dcr; dcc]; 
  ImagePositionPatient = Mvxl2lph(1:3,:)*[1 1 1 1]';

  metadata.PixelSpacing = PixelSpacing;
  metadata.SliceThickness = slcthk;
  metadata.ImageOrientationPatient = ImageOrientationPatient;
  metadata.ImagePositionPatient = ImagePositionPatient + (slice_counter-1)*dcs*slcthk;
  metadata.SliceLocation = (slice_counter-1)*slcthk;

  metadata.InstanceNumber = i;
  metadata.LargestImagePixelValue = maxvol;
  metadata.SmallestImagePixelValue = 0;
  metadata.WindowWidth = (val1-val2);
  metadata.WindowCenter = val3;

  metadata.RescaleSlope = 1;
  metadata.RescaleIntercept = 0;

  if strcmp(metadata.Modality, 'MR')
    metadata.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
    metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
  end

  if strcmp(metadata.Modality, 'CT')
    metadata.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2';
    metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2';
  end

  duid = dicomuid;
  metadata.MediaStorageSOPInstanceUID = duid;
  metadata.SOPInstanceUID = duid;

  fname_out = sprintf('%s/im%3.4i.dcm',output_dir,i);

  switch metadata.BitDepth
    case 8
      metadata.BitsAllocated = 8;
      metadata.BitsStored = 8;
      metadata.HighBit = 7;
      if (multiframe_flag | enhanced_flag) & ~rgb_flag
	dat = reshape(dat, [rows, cols, frames, slices]);
	dicomwrite( uint8( permute(dat, [2 1 3 4]) ), fname_out, metadata, 'CreateMode', 'copy', 'MultiframeSingleFile', 'true' );
      elseif ~rgb_flag
	dicomwrite( uint8( squeeze(dat(:,:,slice_counter,vol_counter))' ), fname_out, metadata, 'CreateMode', 'copy' );
      else
	dicomwrite( uint8( permute(squeeze(dat(:,:,slice_counter,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' );
      end
    case 16
      metadata.BitsAllocated = 16;
      metadata.BitsStored = 16;
      metadata.HighBit = 15;
      if (multiframe_flag | enhanced_flag) & ~rgb_flag
	dat = reshape(dat, [rows, cols, frames, slices]);
	dicomwrite( uint16( permute(dat, [2 1 3 4]) ), fname_out, metadata, 'CreateMode', 'copy', 'MultiframeSingleFile', 'true' );
      elseif ~rgb_flag
	dicomwrite( uint16( squeeze(dat(:,:,slice_counter,vol_counter))' ), fname_out, metadata, 'CreateMode', 'copy' );
      else
	dicomwrite( uint16( permute(squeeze(dat(:,:,slice_counter,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' );
      end
    otherwise
      metadata.BitsAllocated = 16;
      metadata.BitsStored = 16;
      metadata.HighBit = 15;
      if (multiframe_flag | enhanced_flag) & ~rgb_flag
	dat = reshape(dat, [rows, cols, frames, slices]);
	dicomwrite( uint16( permute(dat, [2 1 3 4]) ), fname_out, metadata, 'CreateMode', 'copy', 'MultiframeSingleFile', 'true' );
      elseif ~rgb_flag
	dicomwrite( uint16( squeeze(dat(:,:,slice_counter,vol_counter))' ), fname_out, metadata, 'CreateMode', 'copy' );
      else
	dicomwrite( uint16( permute(squeeze(dat(:,:,slice_counter,:)), [2 1 3]) ), fname_out, metadata, 'CreateMode', 'copy' );
      end
  end

  slice_counter = slice_counter + 1;

  fprintf('%s\n', num2str(i));

end
fprintf('\n%s - Finished \n',mfilename);

end
