clear
dirRoot = 'data0/';
% gr      = gr_readGrid([]);
% gr      = gr_readGrid(gr,[dirRoot],'dr')
% [gr]    = gr_tri(gr);
load data0/gr

usg.h=eb_readHeader([dirRoot '1_usg_.64']);
[usg.h, usg.ts, usg.data] = eb_readTimeStep(usg.h, [1]);

snd.h=eb_readHeader([dirRoot '1_salt.63']);
[snd.h, snd.ts, snd.data] = eb_readTimeStep(snd.h, [1]);

eta.h=eb_readHeader([dirRoot '1_eta2.61']);
[eta.h, eta.ts, eta.data] = eb_readTimeStep(eta.h, [1]);

pts01 = ob_loadPts_01


[ob_eta] = ob_ini(eta.h,pts01);
[ob_eta] = ob_comp_wxy(ob_eta);	
[ob_eta] = ob_comp_wz(ob_eta);
% [ob_eta] = ob_comp_pxy(ob_eta);

[ob_snd] = ob_ini(snd.h,pts01);
[ob_snd] = ob_comp_wxy(ob_snd);	
[ob_snd] = ob_comp_wz(ob_snd);
% [ob_snd] = ob_comp_pxy(ob_snd);

[ob_usg] = ob_ini(usg.h,pts01);
[ob_usg] = ob_comp_wxy(ob_usg);	
[ob_usg] = ob_comp_wz(ob_usg);
% [ob_usg] = ob_comp_pxy(ob_usg);


[obd_usg]=ob_eb2ob_xy(ob_usg,usg)
[obd_snd]=ob_eb2ob_xy(ob_snd,snd)

dd=snd;

