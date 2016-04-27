function [ob_data,t]=ob_observe(ob)
% [ob_data,t]=ob_observe(ob)
% implements generic observation strategy ob for elcirc binary file ob.eb
% ob - observation structure as returned by ob=ob_pts2ob(pts,eb_fname)
% ob_data -- observed data as [nZlevels,nStations,nVect,nTime]
% t - time
%
% Sergey Frolov


% get the data
if strcmp(ob.pts.tType, 'evry')
    tIdx    = [1:ob.eb.h.nSteps];
elseif strcmp(ob.pts.tType, 'last')
    tIdx    = [ob.eb.h.nSteps];
elseif strcmp(ob.pts.tType, 'user')
    tIdx    = ob.pts.itNums;
else
    error(['unknown ob.pts.tType ' sprintf('%s',ob.pts.tType)])
end
[ob.eb.h, ebts, ebd] = eb_readTimeStep(ob.eb.h,tIdx);

% extract for each time step
for i =1:length(tIdx)
    % the if statement makes sure the data is layied out correctly in ob_data
    if ob.eb.h.flagDm == 2
        if ob.eb.h.flagSv == 1
            [ob_data(1,:,1,i)]=ob_eb2ob_xy(ob,ebts{i},ebd(:,:,i));
        elseif ob.eb.h.flagSv == 2
            [ob_data(1,:,1:2,i)]=ob_eb2ob_xy(ob,ebts{i},ebd(:,:,i));
        end
    elseif  ob.eb.h.flagDm == 3
        [ob_data(:,:,:,i)]=ob_eb2ob_xy(ob,ebts{i},ebd(:,:,i));
    end
    t(i,1)  = ebts{i}.t;
end
% again will be cool to check values flagDm flagSv

if strcmp(ob.pts.zType,'var')
	[obd]=ob_eb2ob_z(ob_data,ob);
end