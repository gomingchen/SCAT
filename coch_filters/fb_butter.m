function [B,A] = fb_butter(N,Wn)
% FB_BUTTER  creates a butterworth filterbank using specified input parameters
%
% The MxN output matrix contains M row vectors of time series data N samples long

for i = 1:size(Wn,1)
    [B(i,:),A(i,:)] = butter(N, Wn(i,:));       %% need to address very small bandwidths...  filters do not converge due to roundoff errors, poles/zeros are extremely close to bandpass region
end
