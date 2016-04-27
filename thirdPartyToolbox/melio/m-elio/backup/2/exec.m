function []=exec(cmd,echoflg)
% []=exec(cmd,echoflg)
% executes any executable cmd
% in case cmd fails cmd output is printed on the screen and error is raised
% if echoflg =1 output of the program is prented on the screen
%
% sergey frolov March 11, 2005

if nargin ==1
  [s,w]=system(cmd);
elseif (nargin==2)&(echoflg==1)
  [s,w]=system(cmd,'-echo');
end
if s~=0
   disp(w)
   error(['execution of ' cmd ' failed'])
end

