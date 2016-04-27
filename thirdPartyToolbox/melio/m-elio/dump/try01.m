clear

gr=readGrid([]);
gr=readGrid(gr,'data2','dr')

h_el=readEbHeader('data2/1_elev.61');
[el.h,el.ts,el.data] = readEbTimeStep(h_el,h_el.nSteps);

h_s=readEbHeader('data2/1_salt.63');
[s.h,s.ts,s.data] = readEbTimeStep(h_s,h_s.nSteps);

h_t=readEbHeader('data2/1_salt.63');
[t.h,t.ts,t.data] = readEbTimeStep(h_t,h_t.nSteps);

h_uv=readEbHeader('data2/1_hvel.64');
[uv.h,uv.ts,uv.data] = readEbTimeStep(h_uv,h_uv.nSteps);

h_ww=readEbHeader('data2/1_vert.63');
[ww.h,ww.ts,ww.data] = readEbTimeStep(h_ww,h_ww.nSteps);

hts=read_hotstart('data2/2016_hotstart',gr,1)



for i=0:hts.nlev
    if i == 0
        we(i+1,:) = zeros(1,hts.ne);
    else
        idxl        = find(ww.h.idx.idxLev == i);
        idxn        = ww.h.idx.idxNodes(idxl);
        data_np     = zeros(hts.np,1);
        data_np(idxn) = ww.data(idxl);
        we(i+1,:)   = map_np2ne(gr,data_np)';
    end
end
% hts.we  = we;


dd=ww;
data_np     = zeros(hts.nlev+1,hts.np);
for i =1:hts.nlev
    idxl        = find(dd.h.idx.idxLev == i);
    idxn        = dd.h.idx.idxNodes(idxl);
    data_np(i+1,idxn) = dd.data(idxl,1)';
end
a= data_np(1:end,:);
b= hts.ww(1:end,:);

dd=ssd;
data_np     = zeros(hts.nlev+1,hts.np);
data_sc     = zeros(hts.nlev+1,hts.ns);
for i =1:hts.nlev
    idxl        = find(dd.h.idx.idxLev == i);
    idxn        = dd.h.idx.idxNodes(idxl);
    data_np(i+1,idxn)   = dd.data(idxl,:)';
    data_sc(i+1,:)      = map_np2sc(gr,data_np(i+1,:));
end

a= data_np(2:end,:);
b= hts.ssd(1:end,:);

plot([a(7,:);b(7,:)]')

a= data_sc(2:end,:);
b= hts.ssd(1:end,:);


dd=uv;
eu = zeros(hts.nlev+1,hts.np);
ev = zeros(hts.nlev+1,hts.np);
for i =1:hts.nlev
    idxl        = find(dd.h.idx.idxLev == i);
    idxn        = dd.h.idx.idxNodes(idxl);
    eu(i+1,idxn)  = dd.data(idxl,1)';
    ev(i+1,idxn)  = dd.data(idxl,2)'; 
end

% eu= data_uv_gc(:,:,1);
% ev= data_uv_gc(:,:,2);
hu= hts.uu;
hv= hts.vv;


dd=uv;
data_vn     = zeros(hts.nlev+1,hts.ns);
data_vt     = zeros(hts.nlev+1,hts.ns);
data_tmp_np = zeros(hts.np,2);
data_tmp_ns = zeros(hts.ns,2);
for i =1:hts.nlev
    idxl        = find(dd.h.idx.idxLev == i);
    idxn        = dd.h.idx.idxNodes(idxl);
    data_uv_gc          = data_tmp_np;
    data_uv_gc(idxn,1)  = dd.data(idxl,1);
    data_uv_gc(idxn,2)  = dd.data(idxl,2); 
    data_uv_sc          = data_tmp_ns;
    data_uv_sc(:,1)     = map_np2sc(gr,data_uv_gc(:,1));
    data_uv_sc(:,2)     = map_np2sc(gr,data_uv_gc(:,2));
    data_lc             = map_gc2lc(gr,dd.h,data_uv_sc);
    data_vn(i+1,:)      = data_lc(:,1)';
    data_vt(i+1,:)      = data_lc(:,2)';
end



eb_vn   = data_vn(1:end,:);
eb_vt   = data_vt(1:end,:);
hts_vn  = hts.vn2(1:end,:);
hts_vt  = hts.vt2(1:end,:);

[time,iths,eta1,eta2,we,vn2,vt2,tsd,ssd,peta,ibad,uu,vv,ww,nosm,tnd,snd,q2,xl]...
     =read_hotstart_old('data2/2016_hotstart',gr.hgrid.np,gr.hgrid.ne,gr.sideCent.np,gr.vgrid.nLevels,1);
gr=readGrid([]);
gr=readGrid(gr,'data2/hgrid.gr3','hg');
gr=readGrid(gr,'data2/vgrid.in','vg');
gr=readGrid(gr,'data2/sidecenters_connect.bp','sc');