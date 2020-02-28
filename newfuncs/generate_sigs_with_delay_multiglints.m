function ts = generate_sigs_with_delay_multiglints(delay, tg)

% Generate new signals, with windowed chirp to avoid onset of signals
% Created by Chen Ming on 1/30/2019 
% modified on 4/3/2019 to allow input of parameters
% inputs, delay, the delay between pulse and echo, i.e., between bat and
% target
%         tg, the time delay in microsecond between glints of the target

load(fullfile('newsigs', '1H_100k_20k_fs_500k_3ms_WelchWin.mat'))
brc = Y'.*2;
Fs = 500E3;

t1 = 1500; % in number of samples
noise = zeros(t1,1);
% timeL = 20*1E-3; % in second, the length of the resulted signal
% ts.time = 1/Fs:1/Fs:timeL; % time from 1/Fs to 15 ms
timeLms = round(max(delay)*2/340*1E3 + 11);
timeL = timeLms*1E-3;
ts.time = 1/Fs:1/Fs:timeL; % time from 1/Fs to 15 ms
ts.fs = Fs;

% L = timeL*Fs;

ts.fs = Fs;

L = round(timeL*Fs);
data = zeros(L, length(delay));
%ac_sig = zeros(15E-3*Fs,length(delay));

Ng = floor(tg*1E-6*Fs);
    
for i = 1:length(delay)
    N = floor(delay(i)/340*Fs)*2; % number of samples delayed = delay in distance/340*Fs;
    %echo = brc./2;
    echo = brc./3 + [zeros(Ng,1);brc(1:end-Ng)./3];
    echo_corr = echo/max(abs(echo)).*max(abs(brc)); % corrected on 11/15 to magnify the amplitude of the echo
    data_temp = [noise;brc;zeros(N-length(brc),1);echo_corr;brc(end-Ng+1:end)./3];
    
    len = length(data_temp);
    if len<length(ts.time)
        G = length(ts.time) - len;
        data_temp = [data_temp; zeros(G,1)];
    end
    
    data(:,i) = data_temp;
end

ts.delay = delay;
ts.data = data;
end
