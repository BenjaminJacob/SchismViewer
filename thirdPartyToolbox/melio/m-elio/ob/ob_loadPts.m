function pts=ob_loadPts(fname,zType,tType,itNums)
% pts=ob_loadPts(fname,zType)
% loads point info from fname text file + defines vertical and time strucutre of observations
% zType  - is a description of the vertical structure [var|sfc|btm|non]
%   ['as defined in fname.pts'|'surface'|'bottom'|'none (2d file)']
% tType  - type of time observation ['evry'|'last'|'user']
% itNums - if tType == 'user' itNums is a list [1 2 3 ..] of times tep numbers as defined in the eb file that will be read later
%   !!!attention itNUms will be passed directly to an eb_readTs
%
% Sergey Frolov April 2004

[pts.stName,pts.x,pts.y,pts.z]=textread(fname, '%6c %f %f %f');
pts.numSt   = size(pts.x,1);

%zType
if strcmp(zType,'var')|strcmp(zType,'non') % not implemented strcmp(zType,'sfc')|strcmp(zType,'btm')|
	pts.zType    = zType; 
else
    error(['unknown zType' sprintf('%s',zType)])
end
%tType
if strcmp(tType,'evry')|strcmp(tType,'last')|strcmp(tType,'user') 
	pts.tType    = tType; 
else
    error(['unknown tType' sprintf('%s',tType)])
end
%itNums
if strcmp(tType,'user') & nargin<4
    error('if tType == user you should define itNums')
elseif strcmp(tType,'user')
    pts.itNums    = itNums; 
end