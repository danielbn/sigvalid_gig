function [data_fft_abs , f , data_fft] =  myfft(data,Fs,hann_win)
% this function returns the single-sided spectrum for data 
n = length(data);
if exist('hann_win','var')
    if hann_win == 1
        data=data.*hanning(n)*2;
    end
end
N = n ;% 
% N = 2^(nextpow2(n)); % identify a new input length that is the next power of 2 from the original signal length. This will pad the signal with trailing zeros
data_fft = fft(data,N)/n;
data_fft_abs = abs(data_fft);
data_fft_abs = 2*data_fft_abs(1:N/2+1); % to conserve energy we times 2 the amp
f = Fs/2*linspace(0,1,N/2+1);
