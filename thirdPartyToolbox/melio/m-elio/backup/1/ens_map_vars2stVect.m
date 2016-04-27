function stVect = ens_map_vars2stVect(d,stateInfo,kind)
%stVect = ens_map_vars2stVect(d,stateInfo,[kind])
% maps vars (d.*) to normalized state vector stVect
% stateInfo - holds info on indexes to stVect and normalization coefficients
% kind type of output 'single' or 'double' single is default
% Sergey Frolov, August 2004

if nargin<3
    kind = 'single';
end

if strcmp(kind,'single')
    stVect = single(zeros(stateInfo.idx.fin.all,1));
    stVect(eval(stateInfo.idx.str.e),1) = single(d.e./stateInfo.var.e);
    stVect(eval(stateInfo.idx.str.s),1) = single(d.s./stateInfo.var.s);
    stVect(eval(stateInfo.idx.str.t),1) = single(d.t./stateInfo.var.t);
    stVect(eval(stateInfo.idx.str.u),1) = single(d.u./stateInfo.var.u);
    stVect(eval(stateInfo.idx.str.v),1) = single(d.v./stateInfo.var.v);
elseif strcmp(kind,'double')
    stVect = (zeros(stateInfo.idx.fin.all,1));
    stVect(eval(stateInfo.idx.str.e),1) = (d.e./stateInfo.var.e);
    stVect(eval(stateInfo.idx.str.s),1) = (d.s./stateInfo.var.s);
    stVect(eval(stateInfo.idx.str.t),1) = (d.t./stateInfo.var.t);
    stVect(eval(stateInfo.idx.str.u),1) = (d.u./stateInfo.var.u);
    stVect(eval(stateInfo.idx.str.v),1) = (d.v./stateInfo.var.v);
else
    error('confuesd about the kind parameter')
end
    