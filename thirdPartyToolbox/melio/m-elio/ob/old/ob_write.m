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

fnMiddle    = ob.outInfo.fmiddle;
for i=1:ob.pts.numSt
    if ob.eb.h.flagSv == 1
        fnMiddle    = [fnMiddle '_1'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);
        
    elseif ob.eb.h.flagSv == 2
        fnMiddle    = [fnMiddle '_1'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);
        
        fnMiddle    = [fnMiddle '_2'];
        fname       = fullfile(ob.outInfo.outDir,[ob.outInfo.fpref fnMiddle ob.outInfo.fpost]);        
    end
end