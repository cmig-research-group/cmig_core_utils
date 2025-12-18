function [vol, M, dcminfo] = QD_read_dicomvol(fnames)

ref_info = dicominfo(fnames{1});

% Check for enhanced DICOM
enhanced_flag = 0;
if isfield(ref_info, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end

% Remove CT scout images
ct_flag = 0;
modality = ref_info.Modality;
if ~isempty(regexpi(modality, 'CT'))
  ct_flag = 1;
end

if ct_flag
  counter = 1;
  for i = 1:length(fnames)
    imtype = getfield(dicominfo(fnames{i}),'ImageType');
    if isempty(regexpi(imtype, 'LOCALIZER'))
      fnames2{counter} = fnames{i};
      counter = counter + 1;
    end
  end
  fnames = fnames2;    
end

% Load DICOM header(s), and sort separate files by InstanceNumber
% CCC - Put DICOM header for each file into a field called 'fields' to prevent error when 
%       files from the same volume have different numbers of header tags
%       (Encountered in data exported from MIM, where the first file had additional meta data added
%       to the DICOM header that wasn't present in the other files)  
posvec = NaN(length(fnames), 1);
for fi = 1:length(fnames)
  dcminfo(fi).fields = dicominfo(fnames{fi});
  posvec(fi) = dcminfo(fi).fields.InstanceNumber;
end
[sv, si] = sort(posvec);
fnames = fnames(si);
nfiles = length(fnames);
dcminfo = dcminfo(si);
dcminfo_1 = dcminfo(1).fields;

M = read_dicom_M(fnames);
M = M_LPH_TO_RAS*M; % Convert from DICOM LPH to RAS coordinates

if enhanced_flag % Enhanced DICOM -----------------------------------------------------
  fname = fnames{1};
  vol = squeeze(double(dicomread(fname)));
  vol = permute(vol, [2 1 3]);

  frames = size(vol, 3); % Assumes frame = slice of 3D vol
  for i = 1:frames
    item_str = sprintf('Item_%d', i);
    RescaleSlope = dcminfo_1.PerFrameFunctionalGroupsSequence.(item_str).PixelValueTransformationSequence.Item_1.RescaleSlope;
    RescaleIntercept = dcminfo_1.PerFrameFunctionalGroupsSequence.(item_str).PixelValueTransformationSequence.Item_1.RescaleIntercept;
    vol(:,:,i) = RescaleIntercept + RescaleSlope.*vol(:,:,i);
  end

else % Classic DICOM -----------------------------------------------------------------
  dim1 = double(dcminfo_1.Columns);
  dim2 = double(dcminfo_1.Rows);
  dim3 = nfiles;
  
  vol = zeros(dim1,dim2,dim3);
  for n = 1:nfiles
    fname = fnames{n};
    x = double(dicomread(fname));
    hdr = dcminfo(n).fields;
    if isfield(hdr,'RescaleIntercept') && isfield(hdr,'RescaleSlope')  % AMD: Rescale image values, if fields exist
      x = hdr.RescaleIntercept + hdr.RescaleSlope*x;
    end
    vol(:,:,n) = x';
  end

end

end
