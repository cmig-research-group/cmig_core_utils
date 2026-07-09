function ds = RK_write_segSTRUCT(segSTRUCT,base_dicom, output_dir, SeriesDescription, SeriesNumber)
% function ds = RK_write_segSTRUCT(segSTRUCT,base_dicom, output_dir, SeriesDescription, SeriesNumber)

% set defaults for SeriesNumber and SeriesDescription 
if ~exist('SeriesNumber', 'var'); SeriesNumber = 5003; end
if ~exist('SeriesDescription', 'var'); SeriesDescription = 'RT Structure Set'; end

% create output_dir if it doesn't exist 
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% read base_dicom
if iscell(base_dicom)
    file_list = base_dicom;
else
    file_list = recursive_dir(base_dicom);
end
totfiles = length(file_list);
ivec = nan(1,numel(totfiles));
for i = 1:totfiles;
    fprintf('%s -- %s | File %i of %i\r',datestr(now), mfilename, i, totfiles);
    fnametmp = char(file_list(i));
    try
        dcminfo{i} = dicominfo(fnametmp);
        ivec(i) = dcminfo{i}.InstanceNumber;
    catch
        continue
    end
end

% sort dcminfo  by instance number
[~,ix] = sort(ivec); 
dcminfo = dcminfo(ix); 

% start creating metafile
uid.sop = dicomuid; 
uid.series = dicomuid;

ds.Width = '';
ds.Height = '';
ds.BitDepth = '';
ds.ColorType = '';


ds.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3';
ds.MediaStorageSOPInstanceUID = uid.sop;
ds.TransferSyntaxUID = '1.2.840.10008.1.2';
ds.ImplementationClassUID = '1.2.246.352.70.2.1.160.3';
ds.ImplementationVersionName = 'DCIE 2.2';
ds.SpecificCharacterSet = 'ISO_IR 192';
ds.InstanceCreationDate = datestr(now, 'yyyymmdd'); 
ds.InstanceCreationTime = datestr(now, 'HHMMSS'); 
ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3';
ds.SOPInstanceUID = uid.sop;
ds.StudyDate = dcminfo{1}.StudyDate;
if isfield(dcminfo{1}, 'StudyTime') 
  ds.StudyTime = dcminfo{1}.StudyTime;
end
ds.AccessionNumber = '';
ds.Modality = 'RTSTRUCT';
ds.Manufacturer = 'Unspecified';
ds.ReferringPhysicianName = 'Unspecified';

% Add patient / study identifiers
ds.StationName = 'ro-ariadb-v';
ds.StudyDescription = 'RT Contours';
ds.SeriesDescription = SeriesDescription; ;
ds.PhysiciansOfRecord = 'Unspecified';
ds.OperatorsName = 'Unspecified';
ds.ManufacturerModelName = 'CMIG RT';
ds.PatientName = dcminfo{1}.PatientName; 
ds.PatientID = dcminfo{1}.PatientID;
if isfield(dcminfo{1}, 'PatientBirthDate')
	ds.PatientBirthDate = dcminfo{1}.PatientBirthDate; 
end
if isfield(dcminfo{1}, 'PatientBirthTime')
	ds.PatientBirthTime = dcminfo{1}.PatientBirthTime; 
end
if isfield(dcminfo{1}, 'PatientSex')
  ds.PatientSex = dcminfo{1}.PatientSex;
end 
ds.DeviceSerialNumber = '000000';
ds.SoftwareVersions = '01.01.01';
ds.StudyInstanceUID = dcminfo{1}.StudyInstanceUID; 
ds.SeriesInstanceUID = uid.series;
ds.StudyID = '123456';
if isfield(dcminfo{1}, 'StudyID')
  ds.StudyID = dcminfo{1}.StudyID;
end 
ds.SeriesNumber = SeriesNumber;
ds.InstanceNumber = 1;
ds.StructureSetLabel = sprintf('OAR_%s', datestr(now, 'yyyymmdd')); 
ds.StructureSetDate = datestr(now, 'yyyymmdd'); 
ds.StructureSetTime =  datestr(now, 'HHMMSS');

if ~isfield(dcminfo{1}, 'StudyDescription')
   ds.StudyDescription = 'Unknown Study';
end

% Add Referenced Frame of Reference 
ds = OAR_RS_AddReferencedFrameOfReference(ds, dcminfo); 

% Add Structure Set ROI
ds = OAR_RS_AddStructureSetROI(ds, dcminfo, segSTRUCT); 

% Add Coding scheme 
ds = OAR_RS_AddCodingScheme(ds);

% Add Contours
ds = OAR_RS_AddContour(ds, dcminfo, segSTRUCT);

% Add RT ROI Observation
ds = OAR_RS_AddRTROIObservation(ds, segSTRUCT); 


ds.ApprovalStatus = 'APPROVED';

% write to DICOM
dicomwrite([], sprintf('%s/RS.dcm', output_dir), ds, 'CreateMode', 'Copy'); 

fprintf('\n'); 
return
