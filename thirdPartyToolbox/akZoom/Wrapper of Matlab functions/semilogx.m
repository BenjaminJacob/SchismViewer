function h = semilogx(varargin)
%% wrapper for the original semilogx function
currentFolder = pwd; % save current folder
if ispc
cd([matlabroot '\toolbox\matlab\graph2d']) %go to matlab folder
else
cd([matlabroot '/toolbox/matlab/graph2d']) %go to matlab folder
end

try
  h = semilogx(varargin{:}); % call original function
  akZoom();
  cd(currentFolder) % go back to current folder
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear h;
end
