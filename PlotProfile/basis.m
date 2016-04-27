function nb=basis(b1,b2,u,v)

%transform basis to columnvectors
if size(b1,2)>size(b2,1)
   b1=b1';
   b2=b2';
end

%transform velovities to row vectors
if size(u,1)>size(u,2)
   u=u';
   v=v';
end



                        %u,v, here on column vectors
basis=[b1  b2];        %b1,b2 here on column vectors    
 nb= basis \ [u ; v]; %basisvector length     
 
  %[ xbasis 1  xbasis 2]   [lbasis 1 ]  =       [ X ]
  %[ ybasis 1  y basis2]   [lbasis   ]  =       [ Y ]
 
%   
%  plot(basis(1,1) * nb(1), basis(2,1) * nb(1),'o')
%  plot(basis(1,2) * nb(2), basis(2,2) * nb(2),'x')
 
 
%output nb row1: velocities component of b1
%          row2: velocities component of b12


%u=[1 1 0 -1 -1 -1];
%v=[0 1 1 1 0 -1];

 %b1=[1; 1];
 %b2=[1; -1];
 
 

 