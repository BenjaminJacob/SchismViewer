function ob=ob_ini(pts,eb_fname)
%ob=ob_ini(pts,eb_fname [,outInfo])
% initialises observation structure ob 
% ob binds information on obs points pts, to an elcirc binary file eb, and output directory outInfo 
% ob_ini also compute the info ob.wxy ob.wz for interpolation 
% INPUT/OUTPUT 
% pts=ob_loadPts(fname,zType,...) - strucutre describing observ points
% eb_fname -- string with a name of the file  
%
% Sergey Frolov April 2004

    %inofrmation on points
ob.pts      = pts;
    %ob.outInfo contains info on output configuration, by default no output
ob.outInfo.outFlag = 0;
    % ob.eb - information on the elcirc binary file that will be observed  
ob.eb.h     = eb_readHeader(eb_fname);
    % ob.gr various grid information pertinent to a later interpolation
ob.gr       = ob_tri(ob);
    % ob.wz ob.wxy interpolation indexces and coeffs
[ob]        = ob_get_w(ob);

function gr = ob_tri(ob)
% triangulates horizontal grid in ob.eb.h

	%compute triangulation simmilar to gr_tri
	hgr     = ob.eb.h.hgrid;
	tri     = delaunay(hgr.x,hgr.y);
	hgr.tri = tri;
	if hgr.ne == 0
        hgr.ne     = size(tri,1);
        hgr.elem   = [[1:hgr.ne]', 3*ones(hgr.ne,1), tri, nan*ones(hgr.ne,1)];
	end
	
	gr.hgrid     = hgr;
	gr.vgrid     = ob.eb.h.vgrid;

function ob = ob_get_w(ob)
% compute weights for interpolation on the grid
    
    [ob]=ob_comp_wxy(ob);   %linear interp
%     [ob]=ob_comp_pxy(ob); %nearestr neighbor interp
    if ob.pts.zType == 'var' & ob.eb.h.flagDm > 2
        [ob]=ob_comp_wz(ob);
    end
    %i should check here for the whole range of ob.pts.zType but it not implemented at the moment
    %in case ob.pts.zType == non i will output all layers