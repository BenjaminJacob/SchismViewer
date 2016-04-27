function [ob itgFlag]= ob_ini_staFromPars(gr, xyz)
%[ob]= ob_ini_staFromPars(gr, xyz)
% loads trasect file (*.bp) into ob datastructure
%INPUT/OUTPUT
% gr 	- grid strucutre with hgrid and vgrid fileds
% xyz	- xy vector [x y], z is not implemented 
% ob	- observational data structure
%
% Sergey Frolov, Jun 2005

%initialise ob strucutre
ob      = ob_ini(gr.hgrid, gr.vgrid);

%assign xy points
ob.xy.x     = xyz(:,1);
ob.xy.y     = xyz(:,2);
ob.xy.flag  = 1;    %points loaded
[ob,itgFlag]        = ob_comp_wxy(ob);
ob.xy.flag  = 2;    %weights computed
[ob]=ob_obxy2h(ob);
ob.xy.flag  = 3;    %matrix H assembled
%display('xy matrix H assembled')

