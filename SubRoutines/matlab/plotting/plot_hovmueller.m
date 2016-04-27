function plot_hovmueller(ll,hovmue,dummy1,dummy2,Dim,nrvt,varname)

% function to visualizes of hovmueller diagrams
% based on transect points extracted via 
% hovmue=ExtractHovmueller(fn,gr,ll,varname,days)

fs=14; % fontsize used in plots

if (Dim==62 || Dim == 61)
    nrvti=1;
else
    nrvti=nrvt;
end


localData{1}=load('fort.18_1');
localData{1}(localData{1}==dummy1)=nan;
localData{1}(localData{1}==dummy2)=nan;
localData{1}=localData{1}(:,[1 2:nrvti:end]);
timeSeriesData{1}=localData{1}(:,2:end);
nt=size(timeSeriesData{1},1);


figure(2)
if Dim==62 || Dim==64
    
    
    
    localData2{1}=load('fort.19_1');
    localData2{1}(localData2{1}==dummy1)=nan;
    localData2{1}(localData2{1}==dummy2)=nan;
    localData2{1}=localData2{1}(:,[1 2:nrvti:end]);
    keyboard
    timeSeriesData{1}=timeSeriesData{1}+1i*localData2{1}(:,2:end);
    
    if difference_plot
        localData2{2}=load('fort.19_2');
        localData2{2}(localData2{2}==dummy1)=nan;
        localData2{2}(localData2{2}==dummy2)=nan;
        localData2{2}=localData2{2}(:,[1 2:nrvti:end]);
        timeSeriesData{2}=timeSeriesData{2}+1i*localData2{2}(:,2:end);
    end
    
    %
    %     figure(3)
    %     clf
    %     LegendCell={};
    %     for i=1:1+difference_plot
    %
    %         if difference_plot
    %             subplot(1,2,i)
    %         end
    %         plot(timeSeriesData{i})
    %         hold on
    %         plot(timeSeriesData{i}(1,:),'o');
    %         plot(timeSeriesData{i}(end,:),'x');
    %         xlabel('u')
    %         ylabel('v')
    %
    %         for ii=size(LegendCell,1)+1:size(timeSeriesData,2)*i
    %             LegendCell{ii}=sprintf('P %i',ii);
    %         end
    %
    %         legend(LegendCell)
    %     end
    
    
    figure(2)
    subplot(8,1,1:7)
    imagesc(cumsum(hovmue.ds)/1000,1:nt,real(timeSeriesData{1}))
    ylabel('time / h','fontsize',fs)
    axis tight
    title(sprintf('%s : u at depth %.2f',varname,hovmue.z),'fontsize',fs)
    colorbar('location','NorthOutside')
    subplot(8,1,8)
    plot(cumsum(hovmue.ds)/1000,hovmue.bottom)
    xlabel('profile length','fontsize',fs)
    ylabel('profile depth','fontsize',fs)
    set(gca,'fontsize',fs)
    set(gca,'Ydir','reverse')
    axis tight
    
    
    
    figure(3)
    subplot(8,1,1:7)
    imagesc(cumsum(hovmue.ds)/1000,1:nt,imag(timeSeriesData{1}))
    ylabel('time / h','fontsize',fs)
    axis tight
    title(sprintf('%s : v at depth %.2f',varname,hovmue.z))
    colorbar('location','NorthOutside')
    subplot(8,1,8)
    plot(cumsum(hovmue.ds)/1000,hovmue.bottom)
    xlabel('profile length','fontsize',fs)
    ylabel('profile depth','fontsize',fs)
    set(gca,'fontsize',fs)
    set(gca,'Ydir','reverse')
    axis tight
    
    figure(4)
    subplot(8,1,1:7)
    imagesc(cumsum(hovmue.ds)/1000,1:nt,abs(timeSeriesData{1}))
    ylabel('time / h','fontsize',fs)
    axis tight
    title(sprintf('%s : absolute at depth %f',varname,hovmue.z),'fontsize',fs)
    colorbar('location','NorthOutside')
    subplot(8,1,8)
    plot(cumsum(hovmue.ds)/1000,hovmue.bottom)
    xlabel('profile length','fontsize',fs)
    ylabel('profile depth','fontsize',fs)
    set(gca,'fontsize',fs)
    set(gca,'Ydir','reverse')
    axis tight
    
    
else
    figure(2)
    subplot(8,1,1:7)
    imagesc(cumsum(hovmue.ds)/1000,1:nt,timeSeriesData{1})
    %xlabel('profile length / km')
    ylabel('time / h','fontsize',fs)
    axis tight
    title(sprintf('%s at depth %f',varname,hovmue.z),'fontsize',fs)
    colorbar('location','NorthOutside')
    set(gca,'xtick',[],'fontsize',fs)
    axis tight
    figure(2)
    subplot(8,1,8)
    plot(cumsum(hovmue.ds)/1000,hovmue.bottom)
    xlabel('profile length')
    ylabel('profile depth')
    set(gca,'fontsize',fs)
    set(gca,'Ydir','reverse')
    axis tight
    
end


%make Data available in Workspace
varname(varname=='-')='_';
eval(['workSpaceExport.hovmue.' varname  '=hovmue' ]);
%assignin('base','result',{localData,localData})
assignin('base','result',workSpaceExport)
