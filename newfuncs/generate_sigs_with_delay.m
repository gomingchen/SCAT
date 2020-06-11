function ts = generate_sigs_with_delay(delay)

% Generate new signals, with windowed chirp to avoid onset of signals
% Created by Chen Ming on 1/30/2019 
% modified on 4/3/2019 to allow input of parameters
% input, delay, distance between transmitter and target, in meters

delay = 2.*delay;
load(fullfile('newsigs', '1H_100k_20k_fs_500k_3ms_WelchWin.mat'))
brc = Y'.*2;
Fs = 500E3;

t1 = 1500; % in millisecond
noise = zeros(t1,1);
%timeL = 20*1E-3; % in second, the length of the resulted signal
timeLms = round(max(delay)/340*1E3 + 11); % 11 here is the space after the echo
timeL = timeLms*1E-3;
ts.time = 1/Fs:1/Fs:timeL; % time from 1/Fs to 15 ms
ts.fs = Fs;

L = round(timeL*Fs);

data = zeros(L, length(delay));
%ac_sig = zeros(15E-3*Fs,length(delay));
for i = 1:length(delay)
    N = floor(delay(i)/340*Fs); % number of samples delayed = delay in distance/340*Fs;
    echo = brc./2;
    
    data_temp = [noise;brc;zeros(N-length(brc),1);echo];
    
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
