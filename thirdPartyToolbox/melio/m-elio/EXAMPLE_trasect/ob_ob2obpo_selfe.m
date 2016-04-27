function [obpo]=ob_ob2obpo_selfe(ob, varName, PO, emp, dbInfo, stateInfo)
%[obpo]=ob_ob2obpo_selfe(ob, varName, PO, emp, dbInfo, stateInfo)
% varName is [e | s| ...]
% applies observation operator to projection operator 
% HPO = H*[emp PO]
% Nnc - normalization constants
% Hxe - like HPO just for elev measurments
% Hd  - depth measurments
% 
% Sergey Frolov, June 2005

obsType	= 'xy';
[Hx, ob_mn] = ob_apply_ob2state2(emp, ob, stateInfo, dbInfo, varName, obsType, 0, 0);
[Hxe, junk] = ob_apply_ob2state2(emp, ob, stateInfo, dbInfo, 'e', 'xy',0, 1);
[Hnc, junk] = ob_apply_ob2state2(stateInfo.std.all, ob, stateInfo, dbInfo, varName, obsType,0, 0);
Hd	    = ob.xy.H*ob.gr.hgrid.depth;

numObs		= size(Hx,1);
obpo.obmn	= ob_mn;
obpo.HPO	= zeros(numObs,size(PO,2)+1);
obpo.Hnc        = zeros(numObs,1);
obpo.Hxe        = zeros(ob_mn(1),size(PO,2)+1);

obpo.HPO(:,1)	= Hx;
obpo.Hnc(:,1)   = Hnc;
obpo.Hxe(:,1)   = Hxe;
obpo.Hd         = Hd;

for i =1:size(PO,2)
  [Hx, ob_mn] = ob_apply_ob2state2(PO(:,i), ob, stateInfo, dbInfo, varName, obsType, 0, 0);
  [Hxe, junk] = ob_apply_ob2state2(PO(:,i), ob, stateInfo, dbInfo, 'e', 'xy',0, 1);
  obpo.HPO(:,i+1) = Hx;
  obpo.Hxe(:,i+1) = Hxe;
end


