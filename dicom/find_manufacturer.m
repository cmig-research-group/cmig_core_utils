function manufacturer = find_manufacturer(dcminfo)

patterns_ge = {'GE HEALTHCARE'};
patterns_philips = {'PHILIPS'};
patterns_siemens = {'SIEMENS'};

manufacturer = '';
if isfield(dcminfo, 'Manufacturer')
  manufacturer = dcminfo.Manufacturer;

else
  fprintf('Manufacturer tag not present in DICOM header, attempting to determine manufacturer\n');

  fields = fieldnames(dcminfo);
  for i = 1:length(fields)

    item = dcminfo.(fields{i});
    if ~strcmp(class(item), 'char') && ~strcmp(class(item), 'string')
      continue
    end
    item = char(item);

    match_ge = any(~cellfun(@isempty, regexpi(item, patterns_ge)));
    match_philips = any(~cellfun(@isempty, regexpi(item, patterns_philips)));
    match_siemens = any(~cellfun(@isempty, regexpi(item, patterns_siemens)));
    if match_ge & ~match_philips & ~match_siemens
      manufacturer = 'ge';
      break
    end
    if match_philips & ~match_ge & ~match_siemens
      manufacturer = 'philips';
      break
    end
    if match_siemens & ~match_ge & ~match_philips
      manufacturer = 'siemens';
      break
    end
  end
end

if isempty(manufacturer)
  fprintf('Could not identify manufacturer\n');
  return
end

if ~isempty(regexpi(manufacturer, 'ge'))
  manufacturer = 'ge';
elseif ~isempty(regexpi(manufacturer, 'siemens'))
  manufacturer = 'siemens';
elseif ~isempty(regexpi(manufacturer, 'philips'))
  manufacturer = 'philips';
end

end
