function M_LPH = read_dicom_M(fnames)

if ~iscell(fnames); fnames = {fnames}; end;
nfiles = length(fnames);

ref_info = dicominfo(fnames{1});

% Check for multiframe file
if nfiles == 1
   if isfield(ref_info, 'GridFrameOffsetVector')
      gfov = ref_info.GridFrameOffsetVector;
   end
end

% Check for enhanced DICOM
enhanced_flag = 0;
if isfield(ref_info, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end

dcminfo_1 = ref_info;
if ~enhanced_flag % Traditional DICOMs ------------------------------------------------------

   PixelSpacing = dcminfo_1.PixelSpacing;
   ImageOrientationPatient = dcminfo_1.ImageOrientationPatient;
   ImagePositionPatient_1 = dcminfo_1.ImagePositionPatient;
   
   % Get image volume resolution in slice dimension from DICOM header
   % Will only be used if it can't be computed from ImagePositionPatient values (e.g., multiframe DICOMs with only one image file)
   % First check for SpacingBetweenSlices tag. If not present, fall back to SliceThickness
   if isfield(dcminfo_1, 'SpacingBetweenSlices')
     SpacingBetweenSlices = dcminfo_1.SpacingBetweenSlices;
   else
     SpacingBetweenSlices = dcminfo_1.SliceThickness;
   end
   
   if nfiles > 1
     ImagePositionPatient = NaN(nfiles, 3);
     ipp_diff = NaN(nfiles, 3);
     ipp_diff_norm = NaN(nfiles, 1);
     for i = 1:nfiles
       dcminfo = dicominfo(fnames{i});
       ImagePositionPatient(i,:) = dcminfo.ImagePositionPatient';
       ipp_diff(i,:) = (dcminfo.ImagePositionPatient - ImagePositionPatient_1)';
       ipp_diff_norm(i) = norm(ipp_diff(i,:));
     end

     [unique_ipp, ia, ic] = unique(ImagePositionPatient, 'stable', 'rows');
     unique_ipp_diff_norm = ipp_diff_norm(ia);
     
     indx_end = find(unique_ipp_diff_norm==max(unique_ipp_diff_norm));
     ImagePositionPatient_end = unique_ipp(indx_end,:)';
     steps_to_end = sum(unique_ipp_diff_norm > 0);
     step_along_slice = (ImagePositionPatient_end - ImagePositionPatient_1) / steps_to_end;
   end   

else % Enhanced DICOMs -----------------------------------------------------------------------

  % Get image volume resolution in slice dimension from DICOM header
  % Will only be used if it can't be computed from ImagePositionPatient values (e.g., multiframe DICOMs with only one image file)
  % First check for SpacingBetweenSlices tag. If not present, fall back to SliceThickness
  % (CCC - Not sure how relevant this is for enhanced DICOMs, and I don't have any example datasets without the SpacingBetweenSlices tag)
  if isfield(dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1,'PixelMeasuresSequence')
    PixelSpacing = dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    if isfield(dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1, 'SpacingBetweenSlices')
      SpacingBetweenSlices = dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SpacingBetweenSlices;
    else
      SpacingBetweenSlices = dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
    end
  else
    PixelSpacing = dcminfo_1.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing;
    if isfield(dcminfo_1.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1, 'SpacingBetweenSlices')
      SpacingBetweenSlices = dcminfo_1.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SpacingBetweenSlices;
    else
      SpacingBetweenSlices = dcminfo_1.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
    end
  end

  ImageOrientationPatient = dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient;
  ImagePositionPatient_1 = dcminfo_1.PerFrameFunctionalGroupsSequence.Item_1.PlanePositionSequence.Item_1.ImagePositionPatient;

  frames = dcminfo_1.NumberOfFrames;
  if frames > 1
    ImagePositionPatient = NaN(frames, 3);
    ipp_diff = NaN(frames, 3);
    ipp_diff_norm = NaN(frames, 1);
    for i = 1:frames
      frame_str = sprintf('Item_%d', i);
      ImagePositionPatient(i,:) = dcminfo_1.PerFrameFunctionalGroupsSequence.(frame_str).PlanePositionSequence.Item_1.ImagePositionPatient';
      ipp_diff(i,:) = ImagePositionPatient(i,:) - ImagePositionPatient_1';
      ipp_diff_norm(i) = norm(ipp_diff(i,:));
    end

    [unique_ipp, ia, ic] = unique(ImagePositionPatient, 'stable', 'rows');
    unique_ipp_diff_norm = ipp_diff_norm(ia);
    
    indx_end = find(unique_ipp_diff_norm==max(unique_ipp_diff_norm));
    ImagePositionPatient_end = unique_ipp(indx_end,:)';
    steps_to_end = sum(unique_ipp_diff_norm > 0);
    step_along_slice = (ImagePositionPatient_end - ImagePositionPatient_1) / steps_to_end;
  end
  
end

row_orientation = ImageOrientationPatient(1:3);
col_orientation = ImageOrientationPatient(4:6);
step_along_row = PixelSpacing(2) * row_orientation;
step_along_col = PixelSpacing(1) * col_orientation;

if ~exist('step_along_slice', 'var') % Only a single file/frame
  slice_orientation = cross(row_orientation, col_orientation);
  if exist('gfov', 'var')
    step_along_slice = (gfov(2)-gfov(1)) * slice_orientation;
  else
    step_along_slice = SpacingBetweenSlices * slice_orientation;
  end
end

M = eye(4);
M(1:3,1:3) = [step_along_row step_along_col step_along_slice];
M(1:3,4) = ImagePositionPatient_1 - M(1:3,:)*[1 1 1 1]'; % Adjust for Matlab 1-based indexing

M_LPH = M;

end
