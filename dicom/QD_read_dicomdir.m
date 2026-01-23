function [ctx, dcminfo] = QD_read_dicomdir(dirname, fmt)

if ~exist('fmt','var')
  fmt = '';
end

if ~isempty(fmt)
  flist = dir(sprintf('%s/%s',dirname,fmt));
else
  flist = dir(sprintf('%s/',dirname));
  flist = flist(3:end);
end
fnames = strcat(sprintf('%s/',dirname),{flist.name}');
nfiles = length(fnames);

[vol, M, dcminfo] = QD_read_dicomvol(fnames);

ctx = mgh2ctx(vol, M);

end
