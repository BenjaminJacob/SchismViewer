function [sideCent] = gr_readSCgrid(sccFname,scFname)
%[grid] = gr_readSCgrid(sccFname)
% read connectivity table for side centers
% grid - existing grid structure as outputed by readGrid
% sccFname - side centers connectivity .bp file [ss_node nodeIdx1 nodeIdx2 junk]
% scFname - side centers .bp file [ss_node x y dps]
%
% Sergey Frolov March 08, 2004


[scc]               = gr_readHGrid(sccFname);    
[sc]                = gr_readHGrid(scFname);    
sideCent.nodes      = sc.nodes(:,1:4);
sideCent.connect    = scc.nodes(:,1:3);
sideCent.np         = scc.np;
sideCent.sccFname   = scc.fname;
sideCent.scFname    = sc.fname;
sideCent.ne         = 0;
sideCent.elem       = [];
sideCent.eofLines   = {};
sideCent.nn         = sc.nodes(:,1);
sideCent.x          = sc.nodes(:,2);
sideCent.y          = sc.nodes(:,3);
sideCent.depth      = sc.nodes(:,4);

sideCent.flag      = 1;

    