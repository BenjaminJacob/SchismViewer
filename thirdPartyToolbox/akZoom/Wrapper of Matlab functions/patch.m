function h = patch(varargin)
%% wrapper for the original patch function
currentFolder = pwd; % save current folder

if ispc
cd([matlabroot '\toolbox\matlab\specgraph']) %go to matlab folder
else
cd([matlabroot '/toolbox/matlab/specgraph']) %go to matlab folder
end

try
  h = patch(varargin{:}); % call original function
  cd(currentFolder) % go back to current folder
    akZoom();
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear h;
end