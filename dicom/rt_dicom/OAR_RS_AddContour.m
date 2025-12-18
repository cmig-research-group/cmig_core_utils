function ds = OAR_RS_AddContour(ds, dcminfo, segSTRUCT)

% get colors for the different ROI's 
ROIcolors = OAR_RS_getROIColors(segSTRUCT); 

% Check for enhanced DICOMs
enhanced_flag = 0;
if isfield(dcminfo{1}, 'PerFrameFunctionalGroupsSequence')
  enhanced_flag = 1;
end

for id1 = 1:numel(segSTRUCT)
	tempVol.master_seg = squeeze(segSTRUCT(id1).seg); 
	
	% set up ROIContourSequence
	roi_contour_sequence(id1).ROIDisplayColor = reshape(ROIcolors(id1,:), [], 1); 
	roi_contour_sequence(id1).ContourSequence = [];
	roi_contour_sequence(id1).ReferencedROINumber = segSTRUCT(id1).number;   

	contour_sequence_counter = 1; 

	for id2 = 1:size(tempVol.master_seg, 3)
		tempSlice.im = squeeze(tempVol.master_seg(:,:,id2)); 
		tempSlice.im = 1.0*(tempSlice.im ~= 0); 
		if ~any(tempSlice.im(:))
			continue; 
		end	
		
		% get contours from bwboundaries 
		tempSlice.bw = bwboundaries(tempSlice.im, 'noholes'); 	

		% generate the M matrix 
		if ~enhanced_flag
		  tempSlice.IOP = dcminfo{id2}.ImageOrientationPatient; 
		  tempSlice.IPP = dcminfo{id2}.ImagePositionPatient;
		  tempSlice.PS = dcminfo{id2}.PixelSpacing; 
		else
		  item_str = sprintf('Item_%d', id2);
		  tempSlice.IOP = dcminfo{1}.PerFrameFunctionalGroupsSequence.(item_str).PlaneOrientationSequence.Item_1.ImageOrientationPatient;
		  tempSlice.IPP = dcminfo{1}.PerFrameFunctionalGroupsSequence.(item_str).PlanePositionSequence.Item_1.ImagePositionPatient;
		  tempSlice.PS = dcminfo{1}.PerFrameFunctionalGroupsSequence.(item_str).PixelMeasuresSequence.Item_1.PixelSpacing;
		end
		tempSlice.M = zeros(4,4); 
		tempSlice.M(1:3,1) = tempSlice.IOP(1:3)*tempSlice.PS(1); 
		tempSlice.M(1:3,2) = tempSlice.IOP(4:6)*tempSlice.PS(2); 
		tempSlice.M(1:3,4) = tempSlice.IPP; 
		tempSlice.M(4,4) = 1;

		% iterate through the contours 
		for id3 = 1:numel(tempSlice.bw)
			tempBW.bw = tempSlice.bw{id3}; 
			
			% 0 index
			tempBW.bw = tempBW.bw - 1; 
		
			% contour image sequence 
			contour_image_sequence(1).ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2';
			if ~enhanced_flag
			  contour_image_sequence(1).ReferencedSOPInstanceUID = dcminfo{id2}.SOPInstanceUID; 
			else
			  contour_image_sequence(1).ReferencedSOPInstanceUID = dcminfo{1}.SOPInstanceUID;
			end
			contour_image_sequence = OAR_RS_sequence(contour_image_sequence); 

			% contour data 
			tempBW.cData = cat(2, tempBW.bw(:,1), tempBW.bw(:,2), zeros(size(tempBW.bw, 1), 1), ones(size(tempBW.bw, 1), 1)); 
			tempBW.M_times_cData = tempSlice.M * tempBW.cData'; 
		
			if sum(round(std(tempBW.M_times_cData, [], 2)./0.1)*0.1 == 0) > 2	
				continue; 
			end
	


			tempBW.M_times_cData_col = reshape(tempBW.M_times_cData(1:3,:),[],1); 

			% contour Sequence 
			roi_contour_sequence(id1).ContourSequence(contour_sequence_counter).ContourImageSequence = contour_image_sequence; 
			roi_contour_sequence(id1).ContourSequence(contour_sequence_counter).ContourGeometricType = 'CLOSED_PLANAR'; 
			roi_contour_sequence(id1).ContourSequence(contour_sequence_counter).NumberOfContourPoints = size(tempBW.bw, 1); 
			roi_contour_sequence(id1).ContourSequence(contour_sequence_counter).ContourData = tempBW.M_times_cData_col;  
			contour_sequence_counter = contour_sequence_counter + 1;
			
			clear tempBW
		end
		clear tempSlice
	end
	if ~isempty(roi_contour_sequence(id1).ContourSequence)
		roi_contour_sequence(id1).ContourSequence = OAR_RS_sequence(roi_contour_sequence(id1).ContourSequence); 
	end
end
roi_contour_sequence = OAR_RS_sequence(roi_contour_sequence); 

ds.ROIContourSequence = roi_contour_sequence; 

return
