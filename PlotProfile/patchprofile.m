function ch=patchprofile(nx,nsigma,X1,Y1,data,varargin)



%tri_pt: trinagle vertices  L1 times 3
%rect_pt: rectangle verives L2 times 4
%X1: vertices x-coordinates
%Y1: vertices y-coordinates


%plot from left to right in peaces to avoid horizonta� lines
if nargin==6
n=varargin{1}
else
n=2;
end

if nargin>=7
    Dmin=varargin{2};
    Dmax=varargin{3};
else
    Dmin=min(data);
    Dmax=max(data);
end


%r=max(c)-min(c);
ncol=200;

 cdata=data;
 if nargin >=7
    cdata( find(cdata < Dmin) )= Dmin(1);
    cdata( find(cdata > Dmax) )= Dmax;   
 end


n=1;
for isigma=1:n:nsigma-n
    ii=[(1:n:nx-n)' (1+n:n:nx)' (1+n:n:nx)'  (1:n:nx-n)'];
    li=length(ii);
    jj=[ones(li,1)*isigma ones(li,1)*isigma ones(li,1)*isigma+n...
        ones(li,1)*isigma+n ];
    faces=ii+(jj-1)*nx;
    p_tri = patch('Faces',faces,'Vertices',[X1 Y1],'FaceColor','interp'...
        ,'FaceVertexCData',cdata,'CDataMapping','scaled','EdgeColor','none'); %triangles
    
end
