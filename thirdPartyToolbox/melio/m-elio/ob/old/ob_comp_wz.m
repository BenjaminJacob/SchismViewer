function [ob]=ob_comp_wz(ob)	
% [ob]=ob_comp_wz(ob)
% computes weights for points ob.pts to be interpolated on the vertical grid ob.gr
% uses linear interpolation between bracating layers
%
% Sergey Frolov April 2004

pz      = ob.pts.z;
pnp     = size(pz,1);
vgr     = ob.gr.vgrid.zLevel;

for i=1:pnp
    ilb     = find(diff(vgr < pz(i)) == -1);   %lower bound index
    iub     = ilb +1;   %upper bound index
    w(i,2)  = (pz(i)-vgr(ilb))/(vgr(iub)-vgr(ilb));
    w(i,1)  = 1 - w(i,2);
    zi(i,:) = [ilb iub];
end

ob.wz.zidx     = zi;
ob.wz.w        = w;

