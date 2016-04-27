




varname='salt';


x=1;
y=1;
z=1;
datadir='blubb'

%extract time series
exProgNr=1;
invokeExtract(gr,varname,x,y,z,datadir,exProgNr)


%extract slab
exProgNr=2;
invokeExtract(gr,varname,x,y,z,datadir,exProgNr)

%extract xyt
exProgNr=3;
invokeExtract(gr,varname,x,y,z,datadir,exProgNr)


%extract xyzt
exProgNr=4;
invokeExtract(gr,varname,x,y,z,datadir,exProgNr)


