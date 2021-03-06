function h = imshow(varargin)
%% wrapper for the original image function
currentFolder = pwd; % save current folder
cd([matlabroot '\toolbox\images\imuitools']) %go to matlab folder

try
  h = imshow(varargin{:}); % call original function
  akZoom(gca);
  cd(currentFolder) % go back to current folder
catch err
  cd(currentFolder) % go back to current folder
  rethrow(err)
end

% suppress output if not needed
if nargout == 0
  clear h;
end