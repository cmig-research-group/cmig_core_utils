function QD_ctx_save_mgh(vol_ctx,fname)
% QD_ctx_save_mgh(vol_ctx,fname)
%
% QD_ctx_save_mgh.m - save ctx volume structure as mgh file
%
% Input
%  vol_ctx - 1-based ctx volume structure (double)
%  fname - output mgh filename
%
% created 1/7/2011 by N White


[vol M] = ctx2mgh(vol_ctx);
QD_save_mgh(vol,fname,M);
