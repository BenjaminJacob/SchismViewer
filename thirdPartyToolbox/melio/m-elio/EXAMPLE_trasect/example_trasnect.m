
%load auxilary info
load fnames
load(fnames.grid,'gr');
load(fnames.dbInfo,'dbInfo')
load(fnames.stInfo,'stateInfo')


%generate an interpolator from the buildpoint file
tr_fname='';
[ob]= ob_ini_fromTrasect(gr, tr_fname);
%compute length along the transect
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2+(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen];

%read timestep for corie database at day 3000
[mbi.yyyy, mbi.ww, mbi.dd, mbi.n] = corie_corie2ywdn(3000);
[d, dry]=ens_readVars(stateInfo, dbInfo, mbi);
xfs_now = ens_map_vars2stVect(d, stateInfo, 1);

%interpolate from state format to trasect
[obfs, ob_mn]  = ob_apply_ob2state2(xfs_now, ob, stateInfo, dbInfo, vn, 'xy',0,1);
obfs_d = ob.xy.H*gr.hgrid.depth;
obfs_z = sz_computeZlevels(obfs_d, obfs_el, gr.vgrid);
obfs = reshape(obfs,obpo.obmn);

%plot
xylims=[];
clims=[];
h=myPlotTr(obfs,trLen,obfs_z,xylims,clims);

