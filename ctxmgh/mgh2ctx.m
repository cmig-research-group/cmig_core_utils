function vols_ctx = mgh2ctx(vols, M)

if ~exist('M','var'), M = eye(4); end

dims = size(vols);
width = dims(1); height=dims(2); depth=dims(3); nvol=prod(dims(4:end)); 

[xsize,ysize,zsize,x_r,x_a,x_s,y_r,y_a,y_s,z_r,z_a,z_s,c_r,c_a,c_s] = mat2mgh(M,width,height,depth);

vols_ctx = struct;
vols_ctx.imgs = double(vols);
vols_ctx.Mvxl2lph = M_RAS_TO_LPH * M;
vols_ctx.dimc = height; % size(vol,2);
vols_ctx.dimr = width; % size(vol,1);
vols_ctx.dimd = depth; % size(vol,3);
vols_ctx.vx = xsize;
vols_ctx.vy = ysize;
vols_ctx.vz = zsize;
vols_ctx.lphcent=[-c_r;-c_a;c_s];
vols_ctx.DirCol=[-x_r;-x_a;x_s];
vols_ctx.DirRow=[-y_r;-y_a;y_s];
vols_ctx.DirDep=[-z_r;-z_a;z_s];

end

