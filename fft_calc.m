
function [f, Y, NFFT]=fft_calc(y,Fs)
T = 1/Fs;                     % Sample time
L = length(y);                     % Length of signal
t = (0:L-1)*T;                % Time vector
k=1;
NFFT = 2^nextpow2(L); % Next power of 2 from length of y  -> zero padding could be reason for strange power of M2 frequencies in SPectrum
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
end