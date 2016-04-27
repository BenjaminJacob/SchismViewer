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
ndays=14;
dt=1/24;
tmod=(starttime:dt:starttime+ndays)';

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

    if isfield(stations{i},'slev')
    
    
    pege=pegel{i};
    
    stat.name=pege.name;
    stat.times=pege.times;
    stat.start_time=pege.times(1);
    stat.dt=pege.dt;
    stat.slev.data=pege.elev/100; %conver centimeter to meter
    stat.longitude.data=pege.longitude;
    stat.latitude.data=pege.latitude;
    
    stations{end+1}=stat;
    
    end
    
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

%set minmum depth to prevent choosing to shallow positions
Dmin=2;
lons(gr.hgrid.depth<2)=inf;
lons2(gr2.hgrid.depth<2)=inf;

for istat=1:length(stations)
    plot(stations{istat}.longitude.data,stations{istat}.latitude.data,'kx')
    
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
stack1=3;

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



stack0=1;
stack1=10;

nt2=(stack1-stack0+1)*24;
np=gr2.hgrid.np;

variable='elev.61';

values2=zeros(np,nt2);


for stack=stack0:stack1
    
  h=sz_readHeader(fullfile(fn2.OutDir,sprintf('%i_%s',stack,variable)));
  
  [vals ts]=sz_readTimeStep(h,1:24);
    
  inds=(1:24)+24*(stack-1);
  values2(:,inds)=vals;
    
    
end





values
plot(values(NN,:)')


subplot(10,1,[1 4])
gr_plot(ll.hgrid,ll.hgrid.depth)
hold on
plot(ll2.hgrid.x(domainbd2),ll2.hgrid.y(domainbd2),'k--')
axis tight

for istat=1:length(stations)
    
    subplot(10,1,[1 4])
    hold on 
    ph1=plot(stations{istat}.longitude.data,stations{istat}.latitude.data,'kx','linewidth',4,'Markersize',20)
    title(stations{istat}.name)
    subplot(10,1,[5 10])
    cla
    elev=stations{istat}.slev.data;
    elev(elev<-998)=nan;
    plot(stations{istat}.times,elev,'k')
    hold on
    plot(tmod(1:nt2),values2(NN2(istat),:)','b--')
    plot(tmod(1:nt),values(NN(istat),:)','r--')
    %title(stations{istat}.name)
    legend('gauge','schism old','schism new(21 vert)','location','EastOutside')
    xlim(tmod([1 nt])+[0 3]')
    datetick('x','keeplimits')
    print('-dpng',[stations{istat}.name '_Validation'])
    %pause
    delete(ph1)
end

