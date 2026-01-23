function [vol, M, mr_parms, volsz] = QD_load_mgh(fname,slices,frames,headeronly,keepsingle)
% [vol, M, mr_parms, volsz] = QD_load_mgh(fname,[slices],[frames],[headeronly],[keepsingle])
%
% Required Input:
%   fname: path of the mgh file
% 
% Optional Input:
%   slices: list of one-based slice numbers to load. All
%     slices are loaded if slices is not specified, or
%     if slices is empty, or if slices(1) <= 0.
%     {default = []}
%   frames: list of one-based frame numbers to load. All
%     frames are loaded if frames is not specified, or
%     if frames is empty, or if frames(1) <= 0.
%     {default = []}
%   headeronly: [0|1] whether to load header info only
%     if 1, vol will be empty
%     {default = 0}
%   keepsingle: [0|1] whether to keep float data as single
%     otherwise convert to double
%     {default = 0}
%
% Output:
%   vol: 4D volume (if surface data, size will be [n,1,1,f])
%   M: 4x4 vox2ras transform such that
%     y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
%     indices are 1-based. 
%     If the input has multiple frames,
%     only the first frame is read.
%   mr_parms: [tr flipangle te ti fov]
%   volsz: size(vol). Helpful when using headeronly as vol is [].
%
% See also: fs_save_mgh
%
% copied from freesurfer: 01/18/06
% Prev Mod:               05/02/12 by Don Hagler
% Last Mod:               03/12/18 by Feng Xue
%

vol = [];
M = [];
mr_parms = [];
volsz = [];

if(nargin < 1 | nargin > 5)
  help(mfilename);
  return;
end

max_tmp_strlen = 100;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

% unzip if it is compressed 
if ((length(fname) >=4 && strcmpi(fname((length(fname)-3):length(fname)), '.MGZ')) | ...
		(length(fname) >=3 && strcmpi(fname((length(fname)-2):length(fname)), '.GZ')))
        tempfname = mmil_tempfname;
        new_fname = sprintf('%s.mgh', tempfname);
        cmd = sprintf('zcat < %s > %s', fname, new_fname);
	[s o e] = jsystem(cmd) ; % AMD: hack make make work more reliably on macOS
%        fprintf(1,'file %s written\n',new_fname);
        if s || ~exist(new_fname,'file')
          disp(cmd)
          disp(s)
          disp(o)
          disp(e)
          keyboard
        end
%	fname = temp_fname ;
	fname = new_fname;
	zipflag = 1;
else
	zipflag = 0;
end

if(exist('slices')~=1) slices = []; end
if(isempty(slices)) slices = 0; end
if(slices(1) <= 0) slices = 0; end

if(exist('frames')~=1) frames = []; end
if(isempty(frames)) frames = 0; end
if(frames(1) <= 0) frames = 0; end

if ~exist('headeronly','var') || isempty(headeronly), headeronly = 0; end
if ~exist('keepsingle','var') || isempty(keepsingle), keepsingle = 0; end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
  error('could not open %s for reading',fname);
end
v       = fread(fid, 1, 'int') ; 
ndim1   = fread(fid, 1, 'int') ; 
ndim2   = fread(fid, 1, 'int') ; 
ndim3   = fread(fid, 1, 'int') ; 
nframes = fread(fid, 1, 'int') ;
type    = fread(fid, 1, 'int') ; 
dof     = fread(fid, 1, 'int') ; 
  
if(slices(1) > 0)
  ind = find(slices > ndim3);
  if(~isempty(ind))
    if zipflag, delete(fname); end
    error('some slices exceed nslices');
  end
end

if(frames(1) > 0)
  ind = find(frames > nframes);
  if(~isempty(ind))
    if zipflag, delete(fname); end
    error('some frames exceed nframes');
  end
end

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ; 
if (ras_good_flag)
  delta  = fread(fid, 3, 'float32') ; 
  Mdc    = fread(fid, 9, 'float32') ; 
  Mdc    = reshape(Mdc,[3 3]);
  Pxyz_c = fread(fid, 3, 'float32') ; 

  D = diag(delta);

  Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; % Should this be kept?

  Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;

  M = [Mdc*D Pxyz_0;  ...
	0 0 0 1];
  M = M*[1 0 0 -1; 0 1 0 -1; 0 0 1 -1; 0 0 0 1]; % Convert from 0-based to 1-based indexing
%  ras_xform = [Mdc Pxyz_c; ...
%	0 0 0 1];
  unused_space_size = unused_space_size - USED_SPACE_SIZE ;
else
  warning('header of %s says RAS not good (whatever that means)',fname);
end

fseek(fid, unused_space_size, 'cof') ;
nv = ndim1 * ndim2 * ndim3 * nframes;  
volsz = [ndim1 ndim2 ndim3 nframes];

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

try
  % Determine number of bytes per voxel
  switch type
   case MRI_FLOAT,
    nbytespervox = 4;
   case MRI_UCHAR,
    nbytespervox = 1;
   case MRI_SHORT,
    nbytespervox = 2;
   case MRI_INT,
    nbytespervox = 4;
  end
catch
  if zipflag, delete(fname); end
  error('header of %s has bad type info',fname);
end;

if(headeronly)
  fseek(fid,nv*nbytespervox,'cof');
  if(~feof(fid))
    [mr_parms count] = fread(fid,4,'float32');
    if(count ~= 4) 
%      fprintf('%s: WARNING: error reading MR parms from %s\n',mfilename,fname);
    end
  end
  fclose(fid);
  if zipflag
    delete(fname);
  end
  return;
end

%------------------ Read in the entire volume ----------------%
if(slices(1) <= 0 & frames(1) <= 0)
  switch type % What aboyt type = 2??
   case MRI_FLOAT,
    if keepsingle
      vol = fread(fid, nv, 'float32=>float32');
    else
      vol = fread(fid, nv, 'float32');
    end;
   case MRI_UCHAR,
    vol = fread(fid, nv, 'uchar');
   case MRI_SHORT,
    vol = fread(fid, nv, 'short');
   case MRI_INT,
    vol = fread(fid, nv, 'int');
   case MRI_LONG,
    vol = fread(fid, nv, 'int16'); % hack to work with lq/lqseg.mgz, which is saved with type MRI_LONG, but is actually int16
  end
  if(~feof(fid))
    [mr_parms count] = fread(fid,4,'float32');
    if(count ~= 4) 
%      fprintf('%s: WARNING: error reading MR parms from %s\n',mfilename,fname);
    end
  end
  fclose(fid) ;

  if isempty(vol)
    keyboard
  end

  if zipflag
    delete(fname);
  end
  
  nread = prod(size(vol));
  if(nread ~= nv)
    error('tried to read %d, actually read %d',nv,nread);
  end
  vol = reshape(vol,[ndim1 ndim2 ndim3 nframes]);

  return;
end

%----- only gets here if a subest of slices/frames are to be loaded ---------%


if(frames(1) <= 0) frames = [1:nframes]; end
if(slices(1) <= 0) slices = [1:ndim3]; end

nvslice = ndim1 * ndim2;
nvvol   = ndim1 * ndim2 * ndim3;
filepos0 = ftell(fid);
if keepsingle
  vol = zeros(ndim1,ndim2,length(slices),length(frames),'single');
else
  vol = zeros(ndim1,ndim2,length(slices),length(frames));
end;
nthframe = 1;
for frame = frames

  nthslice = 1;
  for slice = slices
    filepos = ((frame-1)*nvvol + (slice-1)*nvslice)*nbytespervox + filepos0;
    fseek(fid,filepos,'bof');
    
    switch type
     case MRI_FLOAT,
      if keepsingle
        [tmpslice nread]  = fread(fid, nvslice, 'float32=>float32') ; 
      else
        [tmpslice nread]  = fread(fid, nvslice, 'float32') ; 
      end;
     case MRI_UCHAR,
      [tmpslice nread]  = fread(fid, nvslice, 'uchar') ; 
     case MRI_SHORT,
      [tmpslice nread]  = fread(fid, nvslice, 'short') ; 
     case MRI_INT,
      [tmpslice nread]  = fread(fid, nvslice, 'int') ; 
    end

    if(nread ~= nvslice)
      fclose(fid);
      if zipflag
        delete(fname);
      end
      error('reading slice %d, frame %d, tried to read %d, actually read %d',...
        slice,frame,nvslice,nread);
    end

    vol(:,:,nthslice,nthframe) = reshape(tmpslice,[ndim1 ndim2]);
    nthslice = nthslice + 1;
  end

  nthframe = nthframe + 1;
end

% seek to just beyond the last slice/frame %
filepos = (nframes*nvvol)*nbytespervox + filepos0;
fseek(fid,filepos,'bof');

if(~feof(fid))
  [mr_parms count] = fread(fid,5,'float32');
  if(count < 4) 
%    fprintf('%s: WARNING: error reading MR parms from %s\n',mfilename,fname);
  end
end

fclose(fid) ;
if zipflag, delete(fname); end

return;
