function [ob_data]=ob_eb2ob_xy(ob,ebts,ebd)
%[ob_data]=ob_eb2ob_xy(ob,eb)
% implements observation ob.wxy on elcirc data ebd
% ob - observation structure
% [h, ebts, ebd] = eb_readTimeStep(h,n)
% % !! important ebd should be scalar as in [ebIdx,timeIdx], use squeeze(ebd(:,[1|2],:)) for vector
%
% Sergey Frolov, April 2002

% if ob.eb.h.flagDm == 2
%     [ob_data]=ob_2d(ob,ebd);
% elseif ob.eb.h.flagDm == 3
%     [ob_data]=ob_3d(ob,eb.ts,eb.data);
% else
%     error('uknown number of dimension ob.h.flagDm')
% end


if ob.eb.h.flagDm == 2
    if ob.eb.h.flagSv == 1
        [ob_data]=ob_2d(ob,ebd);
    elseif ob.eb.h.flagSv == 2
        [ob_data(:,1)]=ob_2d(ob,ebd(:,1));
        [ob_data(:,2)]=ob_2d(ob,ebd(:,2));        
    end
elseif ob.eb.h.flagDm == 3
    if ob.eb.h.flagSv == 1
        [ob_data]=ob_3d(ob,ebts,ebd);
    elseif ob.eb.h.flagSv == 2
        [ob_data(:,:,1)]=ob_3d(ob,ebts,ebd(:,1));
        [ob_data(:,:,2)]=ob_3d(ob,ebts,ebd(:,2));
    end
else
    error('uknown number of dimension ob.h.flagDm')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%private functions
function [obd]=ob_2d(ob,ebd)
%observation from the 2d field    
    pidx    = ob.wxy.pidx;
    w       = ob.wxy.w;
    obd     = sum(ebd(pidx).*w,2)';    
    
function [obd]=ob_3d(ob,ts,ebd)
%observation from the 3d field
    pidx    = ob.wxy.pidx;
    w       = ob.wxy.w;
    idxLev  = ob.eb.h.idx.idxLev;
    idxNodes= ob.eb.h.idx.idxNodes;
    nlev    = ob.gr.vgrid.nLevels;
    np      = ob.gr.hgrid.np;
    ebidx   = [1:1:np]'; %np+1 will be an index of a nan value
    obd_tmp = nan*ones(1,np+1);
    sidx    = ts.sidx;    
    for i =1:nlev
        idxl        = find(idxLev == i);
        idxn        = idxNodes(idxl);
        idxn(find(sidx < i)) ... %surface index
                    = np+1;
        obdt        = obd_tmp;
        obdt(idxn)  = ebd(idxl)';
        obd(i,:)    = sum(obdt(ebidx(pidx)).*w,2)';
    end

function [obd]=ob_2d_old(ob,ebd)
%observation from the 2d field    
    pidx    = ob.wxy.pidx;
    w       = ob.wxy.w;
    idx_tmp = nan*ones(1,ob.gr.hgrid.np);
    ebidx   = idx_tmp;
    ebidx(ob.eb.h.idx.idxNodes) ...
            = ob.eb.h.idx.idxNodes;
    obd     = sum(ebd(ebidx(pidx)).*w,2);
    