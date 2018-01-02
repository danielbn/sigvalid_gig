function sigg = hipass_filter(data,Fs,cutoff )

if ~exist('Fs','var')
    Fs = 20833;
end

if ~exist('cutoff','var')
    cutoff = 5;
end
L=length(data);  % series length

t= 0:1/Fs:L/Fs- 1/Fs ;  % time scale
% f = Fs/2*linspace(0,1,L/2+1);  % single-sided positive frequency

% data_fft = fft(data,L)/L; % frequency domain 
sigfft =  fft(data,L);
sigfft(1:ceil(cutoff*t(end))) =0;
sigfft(end - ceil(cutoff*t(end))+2:end) = 0;

% sig = ifft(sigfft)*Fs; % back to time domain 
 sigg = ifft(sigfft);
% plot (t,sig); xlabel('Time')