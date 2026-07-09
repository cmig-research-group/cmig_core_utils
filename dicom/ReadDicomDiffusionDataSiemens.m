function [vol, Mvxl2ras, qmat, bvals, dcminfo] = ReadDicomDiffusionDataSiemens(indir);

encoding_problem = 0;

if iscell(indir)
  file_list = indir;
else
  file_list = recursive_dir(indir);
end

totfiles = length(file_list);
ref_file = file_list{1};
ref_info = dicominfo(ref_file);

% Check for enhanced DICOM
enhanced_flag = 0;
if isfield(ref_info, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end


if enhanced_flag % Enhanced DICOM ------------------------------------------------------

  M_LPH = read_dicom_M(ref_file);
  Mvxl2ras = M_LPH_TO_RAS * M_LPH;

  bvals = NaN(totfiles, 1);
  qmat = NaN(totfiles, 3);

  for i = 1:totfiles
      dcminfo = dicominfo(file_list{i});
      if ~isfield(dcminfo, 'Manufacturer')
	dcminfo.Manufacturer = 'SIEMENS'; % Generate tag for other functions, like gradwarp correction
      end

      bval = dcminfo.PerFrameFunctionalGroupsSequence.Item_1.MRDiffusionSequence.Item_1.DiffusionBValue;
      bvals(i) = bval;

      if isfield(dcminfo.PerFrameFunctionalGroupsSequence.Item_1.MRDiffusionSequence.Item_1, 'DiffusionGradientDirectionSequence')
	qvec = dcminfo.PerFrameFunctionalGroupsSequence.Item_1.MRDiffusionSequence.Item_1.DiffusionGradientDirectionSequence.Item_1.DiffusionGradientOrientation';
	qmat(i,:) = qvec;
      end
  end

  rows = ref_info.Columns;
  cols = ref_info.Rows;
  slices = ref_info.NumberOfFrames;
  vol = zeros(rows, cols, slices, totfiles);
  for i = 1:totfiles
    vol_b = double(squeeze(dicomread(file_list{i})));
    vol_b = permute(vol_b, [2 1 3]);
    vol(:,:,:,i) = vol_b;
  end


else % Classic DICOM -----------------------------------------------------------------------

  fnames = {};
  bvals = [];
  qmat = [];
  instances = [];
  slicelocs = [];

  ii = 1;
  for i = 1:totfiles

    fnametmp = file_list{i};

    try 
      dcminfo = dicominfo(fnametmp);
      if ~isfield(dcminfo, 'Manufacturer')
	dcminfo.Manufacturer = 'SIEMENS'; % Generate tag for other functions, like gradwarp correction
      end

      instances(ii) = dcminfo.InstanceNumber;
      slicelocs(ii) = dcminfo.SliceLocation;

      if isfield(dcminfo, 'Private_0019_100c')
	bval = dcminfo.Private_0019_100c;
      elseif isfield(dcminfo, 'B_value')
	bval = dcminfo.B_value;
      elseif isfield(dcminfo, 'Private_0029_1010')
	csa_hdr = spm_dicom_headers12(fnametmp);
	csa_hdr = csa_hdr{1}.CSAImageHeaderInfo;
	csa_fields = {csa_hdr.name};
	ind_b = find(~cellfun(@isempty,regexpi(csa_fields,'b_value')));
	bval_str = {csa_hdr(ind_b).item.val};
	bval_str = bval_str{1};
	bval = str2double(bval_str);
      elseif isfield(dcminfo, 'SequenceName')
	encoding_problem = 1;
	seq_str = dcminfo.SequenceName;
	[ind1, ind2] = regexp(seq_str, '[0-9]+');
	bval_str = seq_str(ind1:ind2);
	bval = str2double(bval_str);
      end

      if length(bval) > 1
	encoding_problem = 1;
	bval = str2double(char(bval'));
      end
      bvals(ii) = bval;

      if bvals(ii) == 0
	qmat(ii,:) = [0 0 0];
      elseif ~encoding_problem && isfield(dcminfo, 'Private_0019_100e')
	qmat(ii,:) = dcminfo.Private_0019_100e;
      elseif ~encoding_problem && isfield(dcminfo, 'DiffusionGradientDirection')
	qmat(ii,:) = dcminfo.DiffusionGradientDirection;
      elseif ~encoding_problem && isfield(dcminfo, 'Private_0029_1010')
	csa_hdr = spm_dicom_headers12(fnametmp);
	csa_hdr = csa_hdr{1}.CSAImageHeaderInfo;
	csa_fields = {csa_hdr.name};
	ind_q = find(~cellfun(@isempty,regexpi(csa_fields,'DiffusionGradientDirection')));
	if isempty(csa_hdr(ind_q).item)
	  qmat(ii,:) = NaN(1,3);
	else
	  qvec = {csa_hdr(ind_q).item.val};
	  qmat(ii,1) = str2double(qvec{1});
	  qmat(ii,2) = str2double(qvec{2});
	  qmat(ii,3) = str2double(qvec{3});
	end
      elseif encoding_problem
	qmat(ii,:) = NaN(1,3);
      end
      
      fnames{ii} = fnametmp;
      ii = ii + 1;

    catch ME
      error(ME.message);
    end

  end


  nrows = double(dcminfo.Rows);
  ncols = double(dcminfo.Columns);
  ps = double(dcminfo.PixelSpacing);

  if isfield(dcminfo,'Private_0019_100a') % Number of images in mosaic
    mosaic_flag = true;
    nslices = double(dcminfo.Private_0019_100a(1));
    nreps = double(length(fnames));
  else % Not mosaic mode
    mosaic_flag = false;
    nslices = length(unique(slicelocs));
    nreps = sum(slicelocs==slicelocs(1));
  end


  [instances, sortinds] = sort(instances);
  fnames = fnames(sortinds);
  bvals = bvals(sortinds);
  qmat = qmat(sortinds,:);
  slicelocs = slicelocs(sortinds);


  bvals_vol = NaN(nreps, 1);
  qmat_vol = NaN(nreps, 3);
  indx_rep = 1;
  for i = 1:nreps
    bvals_vol(i) = bvals(instances==indx_rep);
    qmat_vol(i,:) = qmat(instances==indx_rep,:);
    if mosaic_flag
      indx_rep = indx_rep + 1;
    else
      indx_rep = indx_rep + nslices;
    end
  end
  bvals = bvals_vol;
  qmat = qmat_vol;


  if mosaic_flag == true
    nblocks = ceil(sqrt(nslices));
  else
    nblocks = 1;
  end
  nrows_slice = nrows/nblocks;
  ncols_slice = ncols/nblocks;

  vol = zeros(nrows_slice, ncols_slice, nslices, nreps);
  if mosaic_flag == true
     
    for ithrep = 1:nreps
      im = dicomread(fnames{ithrep});
      if(isempty(im))
	fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
	vol=[];
	return;
      end
      nthslice = 1;
      for ii = 1:nblocks
	for jj = 1:nblocks
          if nthslice > nslices; break;break;end
          rinds = [nrows_slice*(ii-1)+1 nrows_slice*(ii-1)+nrows_slice];
          cinds = [ncols_slice*(jj-1)+1 ncols_slice*(jj-1)+ncols_slice];
          vol(:,:,nthslice,ithrep) = im(rinds(1):rinds(2),cinds(1):cinds(2));
          nthslice = nthslice + 1;
	end
      end
    end

  else % Not mosaic mode

    indx_rep = 1;
    for acq = 1:nreps
      tmp_fnames = fnames(indx_rep:indx_rep+nslices-1);
      for i = 1:nslices
	im = dicomread(tmp_fnames{i});
	if(isempty(im))
	  fprintf('%s: ERROR: could not load pixel data from %s\n',mfilename,fname);
	  vol=[];
	  return
	end
	vol(:,:,i,acq) = double(im);
      end
      indx_rep = indx_rep + nslices;
    end

  end

  vol = permute(vol, [2 1 3 4]); 


  if mosaic_flag == true
    % use spm_dicom_headers to read siemens private CSA header data...
    %    hdr = spm_dicom_headers(fnames{1}); hdr = hdr{1}; % spm8
    hdr = spm_dicom_headers12(fnames{1}); hdr = hdr{1}; % spm12

    % get slice normal vector
    if regexpi(dcminfo.ManufacturerModelName,'TrioTim')
      snvec = str2num([hdr.CSAImageHeaderInfo(24).item(1:3).val]);
    elseif regexpi(dcminfo.ManufacturerModelName,'Skyra')
      snvec = str2num([hdr.CSAImageHeaderInfo(25).item(1:3).val]);
    elseif regexpi(dcminfo.ManufacturerModelName,'Verio')
      snvec = str2num([hdr.CSAImageHeaderInfo(24).item(1:3).val]);
    elseif regexpi(dcminfo.ManufacturerModelName,'Prisma')
      snvec = str2num([hdr.CSAImageHeaderInfo(25).item(1:3).val]);
    else
      fprintf('%s: ERROR: could not load CSA header info for ManufacturerModelName=%s\n',mfilename,dcminfo.ManufacturerModelName);
      vol=[];
      return;
    end

  else

    im1 = dicominfo(fnames{1});
    im2 = dicominfo(fnames{2});
    slicediff = (im2.ImagePositionPatient - im1.ImagePositionPatient)';
    slicediffnorm = slicediff ./ sqrt(sum(slicediff.^2));
    snvec = slicediffnorm;

  end


  M = eye(4);
  M(1:3,1:2) = [dcminfo.PixelSpacing(1)*dcminfo.ImageOrientationPatient(1:3) dcminfo.PixelSpacing(2)*dcminfo.ImageOrientationPatient(4:6)];
  M(1:3,3) = dcminfo.SpacingBetweenSlices*snvec;

  if mosaic_flag == true
    % derive correct ImagePositionPatient from mosaic dimensions.
    ImagePositionPatient = dcminfo.ImagePositionPatient + M(1:3,1:2)*[(ncols-ncols_slice)/2 (nrows-nrows_slice)/2]';
  else
    iminfo = dicominfo(fnames{1});
    ImagePositionPatient = iminfo.ImagePositionPatient;
  end

  M(1:3,4) = ImagePositionPatient - M(1:3,:)*[1 1 1 1]'; % Adjust for Matlab 1-based indexing
  Mvxl2ras = M_LPH_TO_RAS*M; % Convert from DICOM LPH to RAS coordinates

end


% Set q vector for b=0 volumes to 0, rather than NaN
inds_b0 = (bvals == 0);
qmat(inds_b0,:) = zeros(sum(inds_b0), size(qmat,2));


end
