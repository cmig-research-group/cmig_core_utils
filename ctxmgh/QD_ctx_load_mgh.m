function vol_ctx = QD_ctx_load_mgh(fname);
% vol_ctx = QD_ctx_load_mgh(fname);
%
% QD_ctx_load_mgh.m - loads mgh file into ctx volume structure
%
% Input:
%   fname - 0-based 3d, or multiframe mgh (or mgz) filename
%
% Output:
%   vol_ctx - 1-based ctx volume structure 
%
% created on 1/7/2011 by N White


[vol, M] = QD_load_mgh(fname);
vol_ctx = mgh2ctx(vol,M);
