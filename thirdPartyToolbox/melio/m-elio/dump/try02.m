clear
dirRoot='f:/husers/delcirc/da_01/dtests/';

gr=gr_readGrid([]);
gr=gr_readGrid(gr,[dirRoot 'run/'],'dr')


h_ssd=readEbHeader([dirRoot '1_ssd2.63']);
[ssd.h,ssd.ts,ssd.data] = readEbTimeStep(h_ssd,2);

h_tsd=readEbHeader([dirRoot '1_tsd2.63']);
[tsd.h,tsd.ts,tsd.data] = readEbTimeStep(h_tsd,2);

h_usg=readEbHeader([dirRoot '1_usg_.64']);
[usg.h,usg.ts,usg.data] = readEbTimeStep(h_usg,2);

h_el=readEbHeader([dirRoot '1_elev.61']);
[el.h,el.ts,el.data] = readEbTimeStep(h_el,2);

h_eta=readEbHeader([dirRoot '1_eta2.61']);
[eta.h,eta.ts,eta.data] = readEbTimeStep(h_eta,2);

h_snd=readEbHeader([dirRoot '1_salt.63']);
[snd.h,snd.ts,snd.data] = readEbTimeStep(h_snd,h_snd.nSteps);
%%%%%%%%%%%%
hts_new=hotstart_read_da([dirRoot 'run/6_hotstart'],gr,1)
hts_old=hotstart_read_elcirc([dirRoot 'run2/6_hotstart'],gr,1)

usg.h=eb_readHeader([dirRoot 'run/1_usg_.64']);
[usg.h,usg.ts,usg.data] = eb_readTimeStep(usg.h,1);

usg2.h=eb_readEbHeader([dirRoot 'run2/1_usg_.64']);
[usg2.h,usg2.ts,usg2.data] = readEbTimeStep(usg2.h,1);

usg_new=usg;
usg_old=usg2;

[data_lc]=map_gc2lc(gr,usg.h,usg.data);
dd=usg;
eu = zeros(hts_new.nlev+1,hts_new.np);
ev = zeros(hts_new.nlev+1,hts_new.np);
for i =1:hts_new.nlev
    idxl        = find(dd.h.idx.idxLev == i);
    idxn        = dd.h.idx.idxNodes(idxl);
    eu(i+1,idxn)  = data_lc(idxl,1)';
    ev(i+1,idxn)  = data_lc(idxl,2)'; 
%     eu(i+1,idxn)  = dd.data(idxl,1)';
%     ev(i+1,idxn)  = dd.data(idxl,2)';     
end
hu= hts_old.vn2;
hv= hts_old.vt2;

%%%%%%%%%%%
usg3.h=readEbHeader([dirRoot 'run3/1_usg_.64']);
[usg3.h,usg3.ts,usg3.data] = readEbTimeStep(usg3.h,1);

usg4.h=readEbHeader([dirRoot 'run4/1_usg_.64']);
[usg4.h,usg4.ts,usg4.data] = readEbTimeStep(usg4.h,1);
