function []=fort_truncate(fname,nTruncate)
%[]=fort_truncate(fname,nTruncate)
% truncates matrix file fname to nTruncate columns
%
%Sergey Frolov, January 2005

mn = fort_readHead(fname);
if nTruncate > mn(2)
  error(['cant truncate file ' fname ' nTruncate=' nTruncate ' number fo columns in the file ' mn(2)]) 
elseif nTruncate == mn(2)
  return
else
  for i =1:nTruncate
    varOut = fort_readColumn(fname,i,'single');
  end
end

