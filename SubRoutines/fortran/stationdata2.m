
function station=stationdata2(file)

%file=ls('*.nc');
ncid=netcdf.open(file,'NC_NOWRITE');


%list variable names
% varnemes=cell(nvars,1);
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i=1:nvars
% varnames{i} = netcdf.inqVar(ncid,i-1);
% end

%dimensions
% dimname=cell(ndims,1);dimlen=zeros(ndims,1);
% for i=1:ndims
% [dimname{i}, dimlen(i)] = netcdf.inqDim(ncid,i-1);
% end
% dimname
% dimlen



%% station name
namestart= find(file=='/',1,'last')+1;
if isempty(namestart)
    namestart=1;
end
nameend=find(file=='.',1,'last')-1;
station.name=file(namestart:nameend);

namestart= find(station.name=='_',1,'last')+1;
station.name=station.name(namestart:end);

%% time
t=1;
varid = netcdf.inqVarID(ncid,'TIME');
reftime = netcdf.getAtt(ncid,varid,'units');
refY=str2double(reftime(find(reftime=='-',1)-(4:-1:1)));
refM=abs(str2double( reftime(max(find(reftime=='-',2))-(2:-1:1)) ));
refD=abs(str2double( reftime(max(find(reftime=='-',3))-(2:-1:1)) ));
t_offset=datenum(refY,refM,refD);
station.start_time=netcdf.getVar(ncid,varid,t-1,1)+t_offset;%current time
station.times=netcdf.getVar(ncid,varid)+t_offset;%alle times saved in file
station.dt=diff(station.times(1:2));                 %temporal discriti

%plot(station.times,1)
%datetick('x','dd-mmm-yyyy')
%%


varIDs=netcdf.inqVarIDs(ncid);
for varID=varIDs
    varname=lower(netcdf.inqVar(ncid,varID));
    fullname=netcdf.getAtt(ncid,varID,'long_name');
%     if strcmp(varname,'time')
%         continue
%     elseif strmatch('latitude',lower(varname))
%         varname=latitude;
%         elseif strmatch('longitude',lower(varname))
%         varname=longitude;
%     end
%    varname(varname==' ')=[];
 
    eval(sprintf('station.%s.data=netcdf.getVar(ncid,varID)',varname))
    eval(sprintf('station.%s.fullname=fullname',varname));
end
% 
% 
% eval('station.blub=1')
% 
% 
% netcdf.inqAtt(n)
% 
% station.lat=netcdf.getVar(ncid,varid);
% 
% 
% 
% netcdf.inqVar(ncid)
% 
% 
% ID(ncid,'LATITUDE')
% 
% 
% 
% 
% %lon
% varid=netcdf.inqVarID(ncid,'LONGITUDE');
% station.lon=netcdf.getVar(ncid,varid);
% 
% %depth
% varid=netcdf.inqVarID(ncid,'DEPTH');
% station.depth=netcdf.getVar(ncid,varid);
% 
% %temp
% varid=netcdf.inqVarID(ncid,'TEMP');
% station.temp=netcdf.getVar(ncid,varid);
% 
% 
% %elevation
% varid=netcdf.inqVarID(ncid,'SLEV');
% station.elev=netcdf.getVar(ncid,varid);

netcdf.close(ncid)