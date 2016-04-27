clear
% dirRoot='f:/husers/delcirc/da_01/dtests/run/';
dirRoot='data0/';

gr=gr_readGrid([]);
gr=gr_readGrid(gr,[dirRoot],'dr')

usg.h=eb_readHeader([dirRoot '1_usg_.64']);
[usg.h, usg.ts, usg.data] = eb_readTimeStep(usg.h, usg.h.nSteps);

eta.h=eb_readHeader([dirRoot '1_eta2.61']);
[eta.h, eta.ts, eta.data] = eb_readTimeStep(eta.h, eta.h.nSteps);

snd.h=eb_readHeader([dirRoot '1_salt.63']);
[snd.h, snd.ts, snd.data] = eb_readTimeStep(snd.h, snd.h.nSteps);

tnd.h=eb_readHeader([dirRoot '1_temp.63']);
[tnd.h, tnd.ts, tnd.data] = eb_readTimeStep(tnd.h, tnd.h.nSteps);

ssd.h=eb_readHeader([dirRoot '1_ssd2.63']);
[ssd.h, ssd.ts, ssd.data] = eb_readTimeStep(ssd.h, ssd.h.nSteps);

tsd.h = eb_readHeader([dirRoot '1_tsd2.63']);
[tsd.h, tsd.ts, tsd.data] = eb_readTimeStep(tsd.h, tsd.h.nSteps);

kine.h = eb_readHeader([dirRoot '1_kine.63']);
[kine.h, kine.ts, kine.data] = eb_readTimeStep(kine.h, kine.h.nSteps);

mixl.h = eb_readHeader([dirRoot '1_mixl.63']);
[mixl.h, mixl.ts, mixl.data] = eb_readTimeStep(mixl.h, mixl.h.nSteps);

vert.h = eb_readHeader([dirRoot '1_vert.63']);
[vert.h, vert.ts, vert.data] = eb_readTimeStep(vert.h,vert.h.nSteps);

hts_ref=hotstart_read_da([dirRoot '6_hotstart'],gr,1)

%%%%%%%%%%%%


[hts]=mapFile_eb2hts...
    (gr, 'test.hts', 1, eta.ts.t, eta.ts.tit, 1, eta, snd, ssd, tnd, tsd, usg,hts_ref.tnd0,hts_ref.snd0)