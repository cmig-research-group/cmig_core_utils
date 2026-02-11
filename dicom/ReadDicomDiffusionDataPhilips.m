function [vol, Mvxl2ras, qmat, bvals, gwinfo, dcminfo, inds_synth] = ReadDicomDiffusionDataPhilips(indir)

file_list = recursive_dir(indir);
ref_file = file_list{1};
ref_info = dicominfo(ref_file);

[gwinfo, errmsg] = mmil_get_gradwarpinfo(ref_info);

enhanced_flag = 0;
if isfield(ref_info, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end


if enhanced_flag % Enhanced DICOM ------------------------------------------------------

  M_LPH = read_dicom_M(ref_file);
  Mvxl2ras = M_LPH_TO_RAS * M_LPH;
   
  dcminfo = ref_info;

  frames = dcminfo.NumberOfFrames;
  posmat = NaN(frames, 9);
  bvals = NaN(frames, 1);
  qmat = NaN(frames, 3);

  for i = 1:frames
    item_str = sprintf('Item_%d', i);

    ImagePositionPatient = dcminfo.PerFrameFunctionalGroupsSequence.(item_str).PlanePositionSequence.Item_1.ImagePositionPatient;
    ImageOrientationPatient = dcminfo.PerFrameFunctionalGroupsSequence.(item_str).PlaneOrientationSequence.Item_1.ImageOrientationPatient;
    posmat(i,:) = cat(2, ImagePositionPatient', ImageOrientationPatient');

    bval = dcminfo.PerFrameFunctionalGroupsSequence.(item_str).MRDiffusionSequence.Item_1.DiffusionBValue;
    bvals(i) = bval;

    if isfield(dcminfo.PerFrameFunctionalGroupsSequence.(item_str).MRDiffusionSequence.Item_1, 'DiffusionGradientDirectionSequence')
      qvec = dcminfo.PerFrameFunctionalGroupsSequence.(item_str).MRDiffusionSequence.Item_1.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation';
      qmat(i,:) = qvec;
    end

  end

  [C, IA, IC] = unique(posmat,'stable','rows');
  slices = length(C);
  frame_vols = frames / slices;

  vol_unsorted = double(dicomread(ref_file));
  vol_unsorted = permute(vol_unsorted, [2 1 3 4]);
  [rows, cols, ~] = size(vol_unsorted);
  vol = zeros(rows, cols, slices, frame_vols);
  for i = 1:slices
    vol_slice = vol_unsorted(:,:,:,IC==i);
    for j = 1:frame_vols
      vol(:,:,i,j) = vol_slice(:,:,:,j);
    end
  end

  bvals = bvals(IC==1); 
  qmat = qmat(IC==1,:);


else % Classic DICOM -----------------------------------------------------------------------

  totfiles = length(file_list);
  fnames = {};
  instancenumber = [];
  slicelocation = [];
  qmat = [];
  bvals = [];
  iter = 1;
  fprintf('%s -- %s.m:    Reading DICOM headers...\n',datestr(now),mfilename);
  for i = 1:totfiles
    fname = char(file_list(i));

    try
      dcminfo = dicominfo(fname);
      instancenumber(iter) = dcminfo.InstanceNumber;
      slicelocation(iter) = dcminfo.SliceLocation;

      bval = dcminfo.DiffusionBValue;
      bvals = [bvals; bval];
      qvec = dcminfo.DiffusionGradientOrientation';
      qmat = [qmat; qvec];

      fnames{iter} = fname;
      iter = iter+1;

    catch ME
      warning(sprintf('%s - %s\n',ME.message,fname))
      continue
    end

  end

  dcminfo = ref_info; % Sometimes(?) the last file isn't a regular image file, so reset the dicominfo variable

  [instancenumber, sortindx] = sort(instancenumber);
  qmat = qmat(sortindx,:);
  bvals = bvals(sortindx);
  slicelocation = slicelocation(sortindx);
  fnames = fnames(sortindx);

  M_LPH = read_dicom_M(fnames);
  Mvxl2ras = M_LPH_TO_RAS * M_LPH;

  nr = double(dcminfo.Rows);
  nc = double(dcminfo.Columns);
  slices = length(unique(slicelocation));
  nframes = length(instancenumber) / slices;

  fprintf('%s -- %s.m:    Loading DICOM data...\n',datestr(now),mfilename);
  vol = zeros(nr, nc, slices, nframes);
  for i = 1:nframes

    vol_indxs = i:nframes:length(fnames);
    vol_files = fnames(vol_indxs);

    for j = 1:length(vol_files)
      im = double(dicomread(vol_files{j}))';
      finfo = dicominfo(vol_files{j});
      rescale_slope = finfo.RescaleSlope;
      rescale_intercept = finfo.RescaleIntercept;
      scale_slope = double(finfo.MRScaleSlope);
      im = (im * rescale_slope + rescale_intercept) / (rescale_slope * scale_slope);
      vol(:,:,j,i) = im;
    end

  end

  qmat = qmat(1:nframes,:);
  bvals = bvals(1:nframes);

end


% Set q vector for b=0 volumes to 0, rather than NaN
inds_b0 = (bvals == 0);
qmat(inds_b0,:) = zeros(sum(inds_b0), size(qmat,2));

% Check for synthesized volumes
inds_b_nonzero = bvals > 0;
inds_nan = isnan(qmat);
inds_nan = sum(inds_nan, 2) == 3;
inds_q0 = qmat == 0;
inds_q0 = sum(inds_q0, 2) == 3;
inds_synth = inds_b_nonzero & (inds_nan | inds_q0);

if any(inds_synth)
   fprintf('WARNING: DWI data contains synthesized volume(s)\n');
end


end
