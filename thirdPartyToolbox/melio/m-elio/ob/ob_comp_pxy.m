function [ob]=ob_comp_pxy(ob)	
% [ob]=ob_comp_pxy(ob)
% computes idx weights for points ob.pts to be interpolated on grid ob.gr
% ueses nearest point interpolation (not linear interp!!!)
%
% Sergey Frolov April 2004

pxy     = ob.pts.xy;
hgr     = ob.gr.hgrid;
pidx    = dsearch(hgr.x,hgr.y,hgr.tri,pxy(:,1),pxy(:,2));

ob.wxy.pidx     = pidx;
ob.wxy.w        = ones(size(pidx));

