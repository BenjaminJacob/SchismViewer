function [ob]=ob_comp_wxy(ob)	
% [ob]=ob_comp_wxy(ob)
% computes weights for points ob.pts to be interpolated on grid ob.gr
%
% Sergey Frolov April 2004

px      = ob.pts.x;
py      = ob.pts.y;
pnp     = size(px,1);
hgr     = ob.gr.hgrid;
tidx    = tsearch(hgr.x,hgr.y,hgr.tri,px,py);
if find(isnan(tidx))
    error('cant proceed, few of your points are outside of the grid')
end
pidx    = hgr.tri(tidx,:);
%!!! important order points conterclockwise !!!
pidx    = pidx(:,[3 1 2]);  

tx  = hgr.x(pidx);
ty  = hgr.y(pidx);
if pnp==1  %rotate dimensions
    tx = tx';
    ty = ty';
end

%this section of code is borrowed from ELIO v1.33
aum = tx(:,2).* ty(:,3) - tx(:,3).* ty(:,2);
bum = ty(:,2) - ty(:,3);
cum = tx(:,3) - tx(:,2);
ado = tx(:,3).* ty(:,1) - tx(:,1).* ty(:,3);
bdo = ty(:,3) - ty(:,1);
cdo = tx(:,1) - tx(:,3);
atr = tx(:,1).* ty(:,2) - tx(:,2).* ty(:,1);
btr = ty(:,1) - ty(:,2);
ctr = tx(:,2) - tx(:,1);
arei = (aum + ado + atr);
w(:,3) = (atr + btr.* px + ctr.* py)./arei;
w(:,2) = (ado + bdo.* px + cdo.* py)./arei;
w(:,1) = 1.0 - w(:,2) - w(:,3);

ob.wxy.tidx     = tidx;
ob.wxy.pidx     = pidx;
ob.wxy.w        = w;
