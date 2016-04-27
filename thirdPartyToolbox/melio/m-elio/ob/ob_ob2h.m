function H=ob_ob2h(ob)
%H=ob_ob2h(ob)
% converts observation structure ob into an H observation matrix
% H is sparce
% !!!! hadles only 2d varaibles at the moment
%
% Sergey Frolov April 2004

%make indexces of a sparce matrix H(i,j)
i   = ob.wxy.pidx; 
numObs  = size(i,1);
j   = repmat([1:numObs]',1,3);   %assumes that pidx [numObs, 3], where 3 is 3 points of a triangle
i   = i(:);
j   = j(:);
w   = ob.wxy.w(:); 
%dimensions of matrix H=nxm  where 
%m is dimension of a state to observe
%n is a number of measurements
m   = ob.gr.hgrid.np;   
n   = numObs;

H   = sparse(j,i,w,n,m);