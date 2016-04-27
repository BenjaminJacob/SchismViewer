function [d, dry]=ens_readVars(stateInfo, dbInfo, mbi,dryVal)
%[d, dry]=ens_readVars(stateInfo, dbInfo, mbi,dryVal)
% returns data vectors (d.*) and index of dry nodes (dry) from a simulation databse, aplies reading mask (stateInfo.maskRead)
% stateInfo - info on the state, namely fields stateInfo.maskRead.*
% dbInfo    - info about sim database
% mbi       - member index of (ensInfo.memberIdx{1}.yyyy=2002 ...)
% dryVal    - value to put instead of dry nodes "0" by default
% d         - struct with field for each varaible
% dry       - index of 2d and 3d dry nodes
%
% Sergey Frolov, August 2004

if nargin < 4
    dryVal =0;
end

yyyy    = num2str(mbi.yyyy);
ww      = sprintf('%02d',mbi.ww);
dd      = num2str(mbi.dd);

% tic
fname = fullfile(dbInfo.dbDir,[yyyy '-' ww '-' dbInfo.dbVer '/run'],[dd dbInfo.fnames.e]);
dbInfo.fheaders.e   = eb_changeFname(dbInfo.fheaders.e,fname);
[junk, ts, d.e ]  = eb_readTimeStep(dbInfo.fheaders.e,mbi.n);
            %compute dry indexes for 2 d fieled
dry.idx2D   = uint32(find(ts{1}.sidx==0)); 
            %set dry nodes to dryVal (ussually 0)
d.e(dry.idx2D) = dryVal;
d.e = d.e(eval(stateInfo.maskRead.e));


fname = fullfile(dbInfo.dbDir,[yyyy '-' ww '-' dbInfo.dbVer '/run'],[dd dbInfo.fnames.s]);
dbInfo.fheaders.s   = eb_changeFname(dbInfo.fheaders.s,fname);
[junk, ts, d.s] = eb_readTimeStep(dbInfo.fheaders.s,mbi.n);          
            %compute dry indexes for 3 d fieled
dry.idx3D   = uint32(find(d.s==-99)); 
d.s(dry.idx3D) = dryVal;
d.s = d.s(eval(stateInfo.maskRead.s));

fname = fullfile(dbInfo.dbDir,[yyyy '-' ww '-' dbInfo.dbVer '/run'],[dd dbInfo.fnames.t]);
dbInfo.fheaders.t   = eb_changeFname(dbInfo.fheaders.t,fname);
[junk, ts, d.t] = eb_readTimeStep(dbInfo.fheaders.t,mbi.n);          
d.t(dry.idx3D) = dryVal;
d.t = d.t(eval(stateInfo.maskRead.t));

fname = fullfile(dbInfo.dbDir,[yyyy '-' ww '-' dbInfo.dbVer '/run'],[dd dbInfo.fnames.uv]);
dbInfo.fheaders.uv   = eb_changeFname(dbInfo.fheaders.uv,fname);
[junk, ts, duv] = eb_readTimeStep(dbInfo.fheaders.uv,mbi.n);          
duv(dry.idx3D,:) = dryVal;
d.u = duv(eval(stateInfo.maskRead.uv),1);
d.v = duv(eval(stateInfo.maskRead.uv),2);

% toc