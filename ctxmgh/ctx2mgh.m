function [vol, M_RAS] = ctx2mgh(vol_ctx)

vol = single(vol_ctx.imgs);
M_RAS = M_LPH_TO_RAS * vol_ctx.Mvxl2lph;

end
