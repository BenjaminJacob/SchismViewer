function [workMean,workDenom] = ens_meanPrep(d, workMean,workDenom,stateInfo,dryIdx)
%[workMean,workDenom] = ens_meanPrep(d, workMean,workDenom,stateInfo,dryIdx)
% increments mean and denominator
% dryIdx=ensInfo.dryIdx
%
%Sergey frolov, August 2004


workMean(eval(stateInfo.idx.str.e),1) = workMean(eval(stateInfo.idx.str.e),1)+(d.e./stateInfo.var.e);
workMean(eval(stateInfo.idx.str.s),1) = workMean(eval(stateInfo.idx.str.s),1)+(d.s./stateInfo.var.s);
workMean(eval(stateInfo.idx.str.t),1) = workMean(eval(stateInfo.idx.str.t),1)+(d.t./stateInfo.var.t);
workMean(eval(stateInfo.idx.str.u),1) = workMean(eval(stateInfo.idx.str.u),1)+(d.u./stateInfo.var.u);
workMean(eval(stateInfo.idx.str.v),1) = workMean(eval(stateInfo.idx.str.v),1)+(d.v./stateInfo.var.v);

workDenom   = workDenom + 1;
idx = eval(stateInfo.idx.str.e);
workDenom(idx(dryIdx.idx2D)) = workDenom(idx(dryIdx.idx2D))-1;
idx = eval(stateInfo.idx.str.s);
workDenom(idx(dryIdx.idx3D)) = workDenom(idx(dryIdx.idx3D))-1;
idx = eval(stateInfo.idx.str.t);
workDenom(idx(dryIdx.idx3D)) = workDenom(idx(dryIdx.idx3D))-1;
idx = eval(stateInfo.idx.str.u);
workDenom(idx(dryIdx.idx3D)) = workDenom(idx(dryIdx.idx3D))-1;
idx = eval(stateInfo.idx.str.v);
workDenom(idx(dryIdx.idx3D)) = workDenom(idx(dryIdx.idx3D))-1;
