


y=result.elev{1};
Fs=1;
[f, Y, NFFT]=fft_calc(y,Fs);

Amp=abs(Y(1:NFFT/2+1));
Amp(2:end)=2*Amp(2:end);

figure
plot(1./f,Amp)

xlim([0 100])


pwelch(y,[],[],[],1)

r_t_tide

for Pi=1:size(result.Delev_48minRef{1},2)
    
    y=result.Delev_48minRef{1}(:,Pi);
    time=datenum(2011,7,1)+(1+1/24:1/24:30);
    elev=y;
    
    
    
    [name,freq,tidecon,xout]=r_t_tide(time,elev,...
        'interval',1, ...                     % hourly data
        'start',time(1),...               % start time is datestr(tuk_time(1))
        'error','cboot',...                   % coloured boostrap CI
        'synthesis',1);                       % Use SNR=1 for synthesis.
    
    
    
    
    y=result.Delev_240minRef{1}(:,Pi);
    time=datenum(2011,7,1)+(1+1/24:1/24:30);
    elev=y;
    
    
    
    [name2,freq2,tidecon2,xout2]=r_t_tide(time,elev,...
        'interval',1, ...                     % hourly data
        'start',time(1),...               % start time is datestr(tuk_time(1))
        'error','cboot',...                   % coloured boostrap CI
        'synthesis',1);                       % Use SNR=1 for synthesis
    
    
    
    data=[tidecon(:,1) tidecon2(:,1)];
    
    figure
    bar(data)
    x=1:length(tidecon);
    text(x,tidecon2(:,1),name,'rotation',45)
    title(sprintf('P%i',Pi))
    
end