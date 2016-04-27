function [obd]=ob_eb2ob_z(ob_data,ob)
%[ob_data]=ob_eb2ob_xy(ob,eb)
% implements observation 'ob' on data 'ob_data' in vertical dimnsion
% ob - observation structure
% ob_data - data extracted from eb in horizontal domain by ob_eb2ob_xy()
%
% Sergey Frolov, April 2002

if ob.eb.h>flagDm == 2
    error('can;t extract vertical information from a 2D file')
elseif ob.eb.h.flagDm == 3
    if ob.eb.h.flagSv == 1
        [obd]=ob_z(ob_data,ob);
    elseif ob.eb.h.flagSv == 2
        [obd(1,:,1)]=ob_z(ob_data(:,:,1),ob);
        [obd(1,:,2)]=ob_z(ob_data(:,:,2),ob);
    end
else
    error('uknown number of dimension ob.h.flagDm')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%private functions
function [obdo]=ob_z(obd,ob)
%observation 
    [nlevs, nobs] = size(obd);
    
    zidx    = ob.wz.zidx' + repmat([0:nlevs:(nobs-1)*nlevs],size(ob.wz.zidx',1),1);
    w       = ob.wz.w';
    obdo     = sum(obd(zidx).*w,1);