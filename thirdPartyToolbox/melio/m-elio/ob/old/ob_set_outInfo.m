function [ob]=ob_set_outInfo(ob,outDir,fpref,fmiddle,fpost)
%outInfo=ob_outInfo_ini(ob,outDir,[fpref,fmiddle,fpost])
% sets fields of ob.outInfo structure
% outInfo contains information pertinent to the output of ob_data to disk
% ob - observation strucutre
% outDir  - ouput directory
% fpref   - beggining of filename   : like '1_' in 1_snd__1_obd.dat 
%           if not specified use %i_ from ob.eb.h.fname
% fmiddle - middle of file name     : like 'snd_' in 1_snd__1_obd.dat 
%           if not specified use last 4 leters in ob.eb.h.fname like 'snd_' in 1_snd_.63
% fpost   - end of filename         : like '_obd.dat' in 1_snd__1_obd.dat 
%           if not specified use _obd.dat
%
% Sergey Frolov April 2004

if nargin < 2
    error('not enough inputs')
else
	ob.outInfo.outFlag  = 1;
	ob.outInfo.outDir   = outDir;
end
    %
if nargin < 3
    [pathstr,fname,ext,versn]=fileparts(ob.eb.h.fname);
    ob.outInfo.fmiddle    = [num2str(sscanf(fname,'%i_')),'_'];
else
	ob.outInfo.fmiddle    = fmiddle;
end

if nargin < 3
    [pathstr,fname,ext,versn]=fileparts(ob.eb.h.fname);
    ob.outInfo.fpref    = sprintf('%s',fname(end-3:end));
else
	ob.outInfo.fpref    = fpref;
end

if nargin < 5 
    ob.outInfo.fpost = '_obd.dat';
end