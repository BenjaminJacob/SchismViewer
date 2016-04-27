clear
dr_00 = 'f:/filters_01/enkf01_fe06_v01/ens00_m000/run_forecast/';
dr_01 = 'f:/filters_01/enkf01_fe06_v01/ens01_m000/run_forecast/';
% gr      = gr_readGrid([]);
% gr      = gr_readGrid(gr,[dirRoot],'dr')
% [gr]    = gr_tri(gr);
load data0/gr

pts = ob_loadPts('pts_02.pts','non','evry');
ob_eta_00  = ob_ini(pts,[dr_00 '1_eta2.61'])
ob_usg_00  = ob_ini(pts,[dr_00 '1_usg_.64'])
ob_eta_01  = ob_ini(pts,[dr_01 '1_eta2.61'])
ob_usg_01  = ob_ini(pts,[dr_01 '1_usg_.64'])

%[nZlevels,nStations,nVect,nTime]
[obd_eta_00,t]=ob_observe(ob_eta_00);
[obd_usg_00,t]=ob_observe(ob_usg_00);
[obd_eta_01,t]=ob_observe(ob_eta_01);
[obd_usg_01,t]=ob_observe(ob_usg_01);



H=ob_ob2h(ob_eta_00);


plot(squeeze(obd_usg_00(:,7,1,:))',[1:14],'b')


% [ob_eta]=ob_set_outInfo(ob_eta,'outTest/');
% [ob_ssd]=ob_set_outInfo(ob_ssd,'outTest/');
% [ob_usg]=ob_set_outInfo(ob_usg,'outTest/');
% ob_write(ob_eta, obd_eta, t)
