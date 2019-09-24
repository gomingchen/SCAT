function B = fb_firls(N,Wn)
% FB_FIRLS  creates an FIR least squares filterbank using specified input parameters
%
% The MxN output matrix contains M row vectors of time series data N samples long

% setup digital filter parameters
BW = diff(Wn,1,2);
F = [zeros(size(Wn,1),1) Wn(:,1)-BW Wn Wn(:,2)+BW ones(size(Wn,1),1)];  % Cutoff values for BPF
H = [0 0 1 1 0 0];      % Default values for ideal BPF
W = [1 1 1];          % weights for applying LS algorithm

B=zeros(size(Wn,1),N+1);
for i = 1:size(Wn,1)
    B(i,:) = firls(N,F(i,:),H,W);
end
