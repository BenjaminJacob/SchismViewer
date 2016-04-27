function wet = ens_map_dry2wet(dry,stateInfo)
%wet = ens_map_dry2wet(dry,stateInfo)
% maps index of dry nodes to index of wet nodes, includes maskRead 
% dry - dry Idx as returned by [vars,dryIdx]=ens_readVars
% stateInfo - special data structure
% wet - logical array with ones where is wet and 0 otherwise
%
% Sergey Frolov, August 2004


wet = logical(ones(stateInfo.idx.fin.all,1));

idxTmp = eval(stateInfo.idx.str.e);
wet(idxTmp(intersect(eval(stateInfo.maskRead.e),dry.idx2D))) = false;

idxTmp = eval(stateInfo.idx.str.s);
wet(idxTmp(intersect(eval(stateInfo.maskRead.s),dry.idx3D))) = false;

idxTmp = eval(stateInfo.idx.str.t);
wet(idxTmp(intersect(eval(stateInfo.maskRead.t),dry.idx3D))) = false;

idxTmp = eval(stateInfo.idx.str.u);
wet(idxTmp(intersect(eval(stateInfo.maskRead.uv),dry.idx3D))) = false;

idxTmp = eval(stateInfo.idx.str.v);
wet(idxTmp(intersect(eval(stateInfo.maskRead.uv),dry.idx3D))) = false;

