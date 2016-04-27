%% SCHISM Validating.

%% 1) Load SETUP
% 
% %directories
% fn.runDir='/home/benjamin/Desktop/REMOTE/hero_transfer/hero_home/Setups/GermanBightEstuaries/DBucht_Estuaries_enhncd/';
% fn.OutDir='/home/benjamin/Desktop/TRANSFER/waddensea_cuttoff/combined/';
% %fn.tideAna=

fn.runDir='/home/benjamin/Desktop/REMOTE/hero_transfer/hero_home/Setups/GermanBightEstuaries/DBucht_Estuaries_enhncd/';
fn.OutDir='/home/benjamin/Desktop/SELFE_DATA/GERBIHGT_link/correct_timestep/combined/';



%grids
gr.hgrid=gr_readHGrid(fullfile(fn.runDir,'hgrid.gr3'));
ll.hgrid=gr_readHGrid(fullfile(fn.runDir,'hgrid.ll'));


gr2.hgrid=gr_readHGrid(fullfile(fn2.runDir,'hgrid.gr3'));
ll2.hgrid=gr_readHGrid(fullfile(fn2.runDir,'hgrid.ll'));



%get Domain Boundaries
[openbd lanbd domainbd]=plotBound2(gr.hgrid);

%[openbd2 lanbd2 domainbd2]=plotBound2(gr2.hgrid);

%plot
figure
gr_plot_proj(ll.hgrid,ll.hgrid.depth)
hold on
m_plot(ll.hgrid.x(domainbd),ll.hgrid.y(domainbd),'k','linewidth',2)



%integration period
starttime=datenum(2012,6,1);
ndays=15;
dt=1/24;
tmod=(starttime+dt:dt:starttime+ndays)';

%% Load stations

% Load structs of form

%% read PegelOnline Data
fname='W_O_cm_3950020_2012_06_01_2012_08_31.dat';
fname='W_O_cm_9360010_2012_06_01_2012_08_31.dat';
fname='W_O_cm_9510095_2011_01_01_2013_01_01.dat';

pegelDir='/home/benjamin/Desktop/SELFE_Modellierung/Validation/SSE/PegelOnlineDaten';

files=dir(fullfile(pegelDir,'*.dat'))

for ifile=1:length(files)
    
    fname=files(ifile).name;
    fname=fullfile(pegelDir,fname);
    
    pegel_temp=readPegelOnlineDataFile(fname);
    
    if ifile==1;
        pegel=pegel_temp;
        
    else
        pegel(end+(1:length(pegel_temp)))=pegel_temp;
        
    end
end


% Mark positions
hold on
for i=1:length(pegel)
    m_plot(pegel{i}.longitude,pegel{i}.latitude,'kx','linewidth',4,'Markersize',14)
    m_text(pegel{i}.longitude,pegel{i}.latitude,['  ' pegel{i}.name]);
    
end




%    name: 'Whitby'
%     start_time: 7.346869736111111e+05
%          times: [2852x1 double]
%             dt: 0.010416666744277
%            lat: 54.4799995
%            lon: -0.6200000
%          depth: 0v
%           temp: [1x2852 single]
%           elev: [1x2852 single]
% 


%% read MyOcean Stations



%statpath='/media/benjamin/data/Ozeanographie/SELFE Modellierung/gauge_correlation/stations';
statpath='/home/benjamin/Desktop/SELFE_Modellierung/Validation/SSE/all_myoc_stat/201206';



files=dir(fullfile(statpath, '*.nc'));
nfiles=length(files);

p=zeros(nfiles,2);
idel=[];
for i=1:nfiles
    try  stations{i}=stationdata2(fullfile(statpath, files((i)).name));
        % p(i,:)=[stations(i).lon stations(i).lat];%reference points for model
        
        %delete data not in time span
        if stations{i}.times(1) > tmod1 || stations{i}.times(end) < tmod00
            idel=[idel i];
        end
        
        if std(stations{i}.elev) < 1e-5
             idel=[idel i];
        end
            
        
    catch
        fehler(i)=1;
    end
end



%% transform stationdata to pegel data


% for i=1:length(stations)
% 
%     if isfield(stations{i},'slev')
%     
%     stat=stations{i};
%     pege.name=stat.name;
%     pege.times=stat.times;
%     pege.start_time=stat.times(1);
%     pege.dt=stat.dt;
%     pege.elev=stat.slev.data*100; %conver to centimeter
%     pege.longitude=stat.longitude.data;
%     pege.latitude=stat.latitude.data;
%     
%     pegel{end+1}=pege;
%     
%     end
%     
% end


% transform pegeldata to stationdata


for i=1:length(pegel)

   % if isfield(stations{i},'slev')
    
    
    pege=pegel{i};
    stat.source='PegelOn';
    stat.name=pege.name;
    stat.times=pege.times;
    stat.start_time=pege.times(1);
    stat.dt=pege.dt;
    stat.slev.data=pege.elev/100+pege.offsetNHN; %conver centimeter to meter
    stat.longitude.data=pege.longitude;
    stat.latitude.data=pege.latitude;
    
    stations{end+1}=stat;
    
    %end
    
end




hold on
for i=1:length(pegel)
    m_plot(pegel{i}.longitude,pegel{i}.latitude,'kx','linewidth',4,'Markersize',14)
    m_text(double(pegel{i}.longitude),double(pegel{i}.latitude),['  ' pegel{i}.name]);
    
end





% p=zeros(length(files),2);
% npos=zeros(length(files),1);%nextnode
% inarea=npos;


load_stations

%load all stations in folder

% delete station outside domain
idel=[];


minx=min(ll.hgrid.x);
maxx=max(ll.hgrid.x);
miny=min(ll.hgrid.y);
maxy=max(ll.hgrid.y);
idel=[];
for istat=1:length(stations)
    
    
%     [isin ison]=inpolygon(stations{istat}.longitude.data,stations{istat}.latitude.data,...
%         ll.hgrid.x(domainbd),ll.hgrid.y(domainbd));

    [isin ison]=inpolygon(stations{istat}.longitude.data,stations{istat}.latitude.data,...
        [minx maxx maxx minx minx],[miny miny maxy maxy miny ]);
    
    if (~isin & ~ison) || ~isfield(stations{istat},'slev')
        
        idel=[idel istat];
    end
    
    
end

stations(idel)=[];


figure
plot(ll.hgrid.x(domainbd),ll.hgrid.y(domainbd))
hold on

NN=zeros(length(stations),1);
NN2=NN;

lons=ll.hgrid.x;
lats=ll.hgrid.y;
lons2=ll2.hgrid.x;
lats2=ll2.hgrid.y;



binary_file='elev.61';
ifs=1;
stack1=1;
stack2=15;

XYZ=gr.hgrid.nodes(NN,2:end);
[timevec ts_extract]=initialize_fortran_extract(fn,gr,ll,binary_file,ifs,stack1,stack2,XYZ);


%set minmum depth to prevent choosing to shallow positions
Dmin=3;
lons(gr.hgrid.depth<2)=inf;
lons2(gr2.hgrid.depth<2)=inf;

for istat=1:length(stations)
    plot(stations{istat}.longitude.data,stations{istat}.latitude.data,'kx')
    text(double(stations{istat}.longitude.data),double(stations{istat}.latitude.data),['  ' stations{istat}.name])
    
    %find next Neighbour
    [~,NN(istat)]=min(abs(stations{istat}.longitude.data-lons)+...
    abs(stations{istat}.latitude.data-lats));
% 
%     [~,NN2(istat)]=min(abs(stations{istat}.longitude.data-lons2)+...
%     abs(stations{istat}.latitude.data-lats2));
    
    
end

hold on
plot(lons(NN),lats(NN),'ro')
% 
% hold on
% plot(lons2(NN2),lats2(NN2),'kx')
%%

stack0=1;
stack1=15;

nt=(stack1-stack0+1)*24;
np=gr.hgrid.np;

variable='elev.61';

values=zeros(np,nt);




elevs=[];
for stack=stack0:stack1
    
  h=sz_readHeader(fullfile(fn.OutDir,sprintf('%i_%s',stack,variable)));
  
  [vals ts]=sz_readTimeStep(h,1:24);
    
  inds=(1:24)+24*(stack-1);
  values(:,inds)=vals;
    
    
end


 [xs ys]=ginput(1);
 [was wo]=min(abs(ll.hgrid.x-xs)+abs(ll.hgrid.y-ys));
 %NN for Delfzijl 5424
 sz_readTimeseries_pos2
 

 extract_dir='/home/benjamin/Desktop/SCHISM/Post-Processing-Fortran/outputlink';

 fid=fopen(fullfile(extract_dir,'station.sta'),'w')
 fprintf(fid,'1 1 0 !itransect,ics,ifs (itransect=0: time series; =1: transect. ics=1: Cartesian coordinates; =2: lat/lon. ifs=0: input z are z coordinates; =1: input z are relative to free surface; and in this case, make sure all z>=0\n')
 fprintf(fid,'1 !total # of points\n');
 fprintf(fid,'%f %f 0 !point #, x, y, z',gr.hgrid.x(wo),gr.hgrid.y(wo));
 fclose(fid);
 A=load(fullfile(extract_dir,'fort.18'));
 
 figure
 istat=16;
 elev=stations{istat}.slev.data;
    elev(elev<-998)=nan;
    plot(stations{istat}.times-stations{istat}.utc/24,elev,'k')
    hold on 
    plot(tmod(1:nt),A(:,2),'r--')
    
    %defizil 9086
    %wittduen 247131
    %Viada
    
    B=load('/home/benjamin/Desktop/SCHISM/Post-Processing-Fortran/ts_extract/fort.18');
    
% 
% stack0=1;
% stack1=10;
% 
% nt2=(stack1-stack0+1)*24;
% np=gr2.hgrid.np;
% 
% variable='elev.61';
% 
% values2=zeros(np,nt2);
% 
% 
% for stack=stack0:stack1
%     
%   h=sz_readHeader(fullfile(fn2.OutDir,sprintf('%i_%s',stack,variable)));
%   
%   [vals ts]=sz_readTimeStep(h,1:24);
%     
%   inds=(1:24)+24*(stack-1);
%   values2(:,inds)=vals;
%     
%     
% end
% 
% 
% 


values
plot(values(NN,:)')


subplot(10,1,[1 4])
gr_plot(ll.hgrid,ll.hgrid.depth)
hold on
%plot(ll2.hgrid.x(domainbd2),ll2.hgrid.y(domainbd2),'k--')
axis tight

for istat=1:length(stations)
    
    subplot(10,1,[1 4])
    hold on 
    ph1=plot(stations{istat}.longitude.data,stations{istat}.latitude.data,'kx','linewidth',4,'Markersize',20);
    if isfield(stations{istat},'source');
        stations{istat}.utc=+1;
        title([stations{istat}.name ' (' stations{istat}.source  ')'])
    else
        stations{istat}.utc=0;
        title(stations{istat}.name)
    end
    subplot(10,1,[5 10])
    cla
    elev=stations{istat}.slev.data;
    elev(elev<-998)=nan;
    plot(stations{istat}.times-stations{istat}.utc/24,elev,'k')
    hold on
    %plot(tmod(1:nt2),values2(NN2(istat),:)','b--')
    plot(tmod(1:nt),values(NN(istat),:)','r--')
    %title(stations{istat}.name)
    legend('gauge','schism new(21 vert)','location','EastOutside')
    xlim(tmod([1 nt])+[0 0.2]')
    datetick('x','keeplimits')
    print('-dpng',[stations{istat}.name '_Validation'])
    %pause
    delete(ph1)
end


%%
sim_stations=cell(size(stations));
for i=1:length(stations)

    sim_stations{i}.name=stations{i};
    sim_stations{i}.dt=1;
    sim_stations{i}.longitude.data=ll.hgrid.x(NN(i));
    sim_stations{i}.latitude.data=ll.hgrid.y(NN(i));
    sim_stations{i}.times=tmod;
    sim_stations{i}.slev.data=values(NN(i),:);
    
end

elev=stations{1}.slev.data;
elev(elev<-998)=nan;
plot(stations{1}.times-stations{1}.utc/24,elev,'k')



hold on
elev=stations{24}.slev.data;
elev(elev<-998)=nan;
plot(stations{24}.times-1/24,elev,'g--')



stations{4}.longitude
stations{21}.longitude