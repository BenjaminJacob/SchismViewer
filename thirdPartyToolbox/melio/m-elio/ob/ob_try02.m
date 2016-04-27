clear
dirRoot = 'data0/';
% gr      = gr_readGrid([]);
% gr      = gr_readGrid(gr,[dirRoot],'dr')
% [gr]    = gr_tri(gr);
load data0/gr

pts = ob_loadPts('pts_01.pts','non','evry');
ob_eta  = ob_ini(pts,[dirRoot '1_eta2.61'])
ob_ssd  = ob_ini(pts,[dirRoot '1_ssd2.63'])
ob_usg  = ob_ini(pts,[dirRoot '1_usg_.64'])

[ob_eta]=ob_set_outInfo(ob_eta,'outTest/');
[ob_ssd]=ob_set_outInfo(ob_ssd,'outTest/');
[ob_usg]=ob_set_outInfo(ob_usg,'outTest/');

[obd_eta,t]=ob_observe(ob_eta);
[obd_ssd,t]=ob_observe(ob_ssd);
[obd_usg,t]=ob_observe(ob_usg);

ob_write(ob_eta, obd_eta, t)