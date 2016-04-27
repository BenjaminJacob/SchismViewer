

%% Load Station Data
%stations which lie in domain

statpath='/home/benjamin/Desktop/SELFE_Modellierung/Validation/SSE/all_myoc_stat/201206';

files=dir(fullfile(statpath, '*201206*.nc'));
nfiles=length(files);
p=zeros(length(files),2);
npos=zeros(length(files),1);%nextnode
inarea=npos;

p=zeros(nfiles,2);
for i=1:nfiles
    fname=files((i)).name;
    try
    stations{i}=stationdata(fullfile(statpath, fname));
   % p(i,:)=[stations{i}.lon stations{i}.lat];%reference points for model
   stations{i}.name=fname(find(fname=='_',1,'last')+1:...
       find(fname=='.',1,'last')-1);
    catch
    fehler(i)=1;
    end
end

%% Change Order
% 
% stations_sort=stations_sort_CCW(stations);
% 
% a=stations{6}
% stations{6}=stations{7}
% stations{7}=a
% 
% gr_plot(
% 
% 
% figure
% clf
% 
% m_proj('MERCATOR','long',[-5 10],'lat',[48 60])
% 
% 
% 
% meanx=0;
% meany=0;
% 
% 
% nstations=length(stations)-1
% m_coast
% for i=1:nstations
%     
%     
%     hold on
%     m_plot(double(stations{i}.lon),double(stations{i}.lat),'k.')
%     m_text(double(stations{i}.lon),double(stations{i}.lat),num2str(i))
%     
% %     meanx=meanx+double(stations{i}.lon);
% %     meany=meany+double(stations{i}.lat)
% %     
% %     lons(i)=double(stations{i}.lon);
% %     lats(i)=double(stations{i}.lat);
%     
% end
% 
% leg=[];
% for i=1:nstations
% leg=[leg sprintf('%i: %s ',i,stations{i}.name)    ]
% end
% 1: Whitby 2: Cromer 3: Lowestoft 4: Sheerness 5: Dover 6: HoekVanHolland 7: DenHelder 8: Harlingen 9: Huibertgat 10: Borkum 11: Emden 12: Wilhelmshaven 13: Bremerhaven 14: Cuxhaven 15: Buesum 16: Helgoland 17: Wittduen 18: Vidaa 19: HavnebyRomo 20: Esbjerg 21: HvideSandeYdermole 22: TorsmindeYdermole 23: FerringHavn 24: ThyboronHavn 
% 
% 
% hold on
% m_plot(6.587638889,54.148611111,'ko','Markerfacecolor','y');
% m_text(6.587638889,54.148611111,'Fino1');
% 
% 
% hold on
% m_plot(7.158333333,55.195,'ko','Markerfacecolor','y');
% m_text(7.158333333,55.195,'Fino3');
% m_grid('box','fancy','tickdir','in')
% 
% 
% figure
% gr_plot_proj(gr.hgrid,gr.hgrid.depth,[],[-4 10.5 48 58])
% m_gshhs_f('patch',[.7 .7 .7])
% 
% print(gcf,'-dpng','Overview.png')
% 
