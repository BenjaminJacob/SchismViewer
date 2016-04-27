function h = plot(varargin)
%% wrapper for the original plot function
currentFolder = pwd; % save current folder
if ispc
cd([matlabroot '\toolbox\matlab\graph2d']) %go to matlab folder^M
else
    cd([matlabroot '/toolbox/matlab/graph2d']) %go to matlab folder^M %go to matlab folder
end
try
  h = plot(varargin{:}); % call original function
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
