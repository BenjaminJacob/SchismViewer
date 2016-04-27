function wetIdx = ens_map_dry2wet(dry,stateInfo)
%wet = ens_map_dry2wet(dry,stateInfo)
% maps index of dry nodes to index of wet nodes, includes maskRead 
% dry - dry Idx as returned by [vars,dryIdx]=ens_readVars
% stateInfo - special data structure
% wet - index of wet nodes in a state vector
%!!! this code will be flowed if mask read is non-trivial (some areas ommitied)
%
% Sergey Frolov, August 2004


% wet = logical(ones(stateInfo.idx.fin.all,1));
wetIdx=[];

idxTmp = eval(stateInfo.idx.str.e)';
wetIdx = [wetIdx; idxTmp(setdiff(eval(stateInfo.maskRead.e),dry.idx2D))];

idxTmp = eval(stateInfo.idx.str.s)';
wetIdx = [wetIdx; idxTmp(setdiff(eval(stateInfo.maskRead.s),dry.idx3D))];

idxTmp = eval(stateInfo.idx.str.t)';
wetIdx = [wetIdx; idxTmp(setdiff(eval(stateInfo.maskRead.t),dry.idx3D))];

idxTmp = eval(stateInfo.idx.str.u)';
wetIdx = [wetIdx; idxTmp(setdiff(eval(stateInfo.maskRead.uv),dry.idx3D))];

idxTmp = eval(stateInfo.idx.str.v)';
wetIdx = [wetIdx; idxTmp(setdiff(eval(stateInfo.maskRead.uv),dry.idx3D))];

