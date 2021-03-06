function h = image(varargin)
%% wrapper for the original image function
currentFolder = pwd; % save current folder
if ispc
cd([matlabroot '\toolbox\matlab\specgraph']) %go to matlab folder
else
    cd([matlabroot '/toolbox/matlab/specgraph']) %go to matlab folder
end
try
  h = image(varargin{:}); % call original function
  %akZoom(gca);
  cd(currentFolder) % go back to current folder
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear h;
end