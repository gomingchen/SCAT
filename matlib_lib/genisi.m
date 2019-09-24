function I = genisi(S,fs)
% GENISI  generates the interspike interval (ISI) of the given spike train
%
% I = GENISI(SPIKES,FS) returns a vector of the time differences between
% spikes in SPIKES.
%

spikeLoc = find(S > 0);
spikeTime = spikeLoc/fs;
I = diff(spikeTime);