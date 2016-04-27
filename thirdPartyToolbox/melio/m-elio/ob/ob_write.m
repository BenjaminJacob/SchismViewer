function []=ob_write(ob, obd, t)
% []=ob_write(ob)
% write observation data obd to disk as defined in ob
% ob - observation structure as returned by ob=ob_pts2ob(pts,eb_fname)
% ob_data -- observed data
%
% Sergey Frolov, April 2004

    %do i need to output?
if ~ob.outInfo.outFlag
  display('no input requaeted exiting ob_write')  
  return
end
    % is ob.outInfo.outDir a directory?
if exist(ob.outInfo.outDir)~= 7
    error([sprintf('%s',ob.outInfo.outDir) ' should be a directory'])
end


    % loop through stations and output a file for each one of them
    % all iffs are for correct handling of multidimensional arrays
for i=1:ob.pts.numSt
    if ob.eb.h.flagSv == 1
        fnMiddle    = [ob.pts.stName(i,:) '_' ob.outInfo.fmiddle '_d1'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);
        if ob.eb.h.flagDm == 2
            outVar      = [t squeeze(obd(:,i,1,:))];
        else ob.eb.h.flagDm == 3
            outVar      = [t squeeze(obd(:,i,1,:))'];            
        end
        save(fname,'outVar','-ascii');
    elseif ob.eb.h.flagSv == 2
        fnMiddle    = [ob.pts.stName(i,:) '_' ob.outInfo.fmiddle '_d1'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);
        if ob.eb.h.flagDm == 2
            outVar      = [t squeeze(obd(:,i,1,:))];
        else ob.eb.h.flagDm == 3
            outVar      = [t squeeze(obd(:,i,1,:))'];            
        end
        save(fname,'outVar','-ascii');

        fnMiddle    = [ob.pts.stName(i,:) '_' ob.outInfo.fmiddle '_d2'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);        
        if ob.eb.h.flagDm == 2
            outVar      = [t squeeze(obd(:,i,2,:))];
        else ob.eb.h.flagDm == 3
            outVar      = [t squeeze(obd(:,i,2,:))'];            
        end
        save(fname,'outVar','-ascii');
    end
end