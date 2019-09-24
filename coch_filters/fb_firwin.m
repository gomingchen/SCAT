function B = fb_firwin(N,Wn)
% FB_FIRWIN  creates an FIR filterbank using the frequency windowing method
%
% The MxN output matrix contains M row vectors of time series data N samples long

% setup digital filter parameters
%BW = diff(Wn,1,2);
%F = [zeros(size(Wn,1),1) Wn(:,1)-BW Wn Wn(:,2)+BW ones(size(Wn,1),1)];  % Cutoff values for BPF

B=zeros(size(Wn,1),N+1);
W=chebwin(N+1,80);        % window method to use
for i = 1:size(Wn,1)
    B(i,:) = fir1(N,Wn(i,:),W);
end
